#!/usr/bin/env python3
"""
TeraFlux: A Python library for 
Thermodynamically Enabled and Reaction Attuned FLUXome estimation.

This script integrates gene expression data (FPKM) with a genome-scale metabolic
model (GEM) to predict metabolic fluxes. It uses the principle of maximum entropy
production, formulated as a non-linear optimization problem and solved with CasADi
and IPOPT.
"""

import argparse
import gc
import os
import re
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Tuple, Any

import casadi as ca
import cobra
import numpy as np
import pandas as pd
import psutil

# --- Constants ---
# Define constants for readability and easier modification.
REACTION_DEFAULT_UPPER_BOUND =  ca.inf #30000.0
REACTION_DEFAULT_LOWER_BOUND = -ca.inf #-30000.0
EPSILON = 0  # Small value to prevent log(0) errors in fluxes
EPSILON_g = 1e-6  # Small value to prevent log(0) errors in gene expression

# A pre-compiled regular expression for finding gene IDs in GPR strings.
# This is much faster than repeated string manipulations.
GENE_ID_RE = re.compile(r'\b[a-zA-Z0-9_.-]+\b')


def parse_gpr(rule: str) -> List[List[str]]:
    """
    Parses a Gene-Protein-Reaction (GPR) rule into a list of "AND" groups.
    This structure represents the rule in Disjunctive Normal Form (OR of ANDs).

    Example:
        '(g1 and g2) or g3' is parsed to [['g1', 'g2'], ['g3']]

    Args:
        rule (str): The GPR string from the model.

    Returns:
        List[List[str]]: A list of lists, where each inner list represents a
                         complex (genes connected by "AND") and the outer list
                         represents isozymes (complexes connected by "OR").
    """
    if not rule:
        return []
    # Split the rule by "or" to get the isozyme groups.
    or_groups = []
    for sub_rule in rule.split(" or "):
        # Find all gene IDs within the sub-rule (the "AND" part).
        # This regex handles gene names that may include dots or dashes.
        genes = GENE_ID_RE.findall(sub_rule)
        if genes:
            or_groups.append(genes)
    return or_groups


def preprocess_reaction_data(model: cobra.Model, fpkm_df: pd.DataFrame) -> pd.DataFrame:
    """
    Pre-processes FPKM data and maps it to reactions in the model.

    It performs the expensive GPR parsing
    and 'g' value calculation only ONCE per model, instead of on every iteration.

    Args:
        model (cobra.Model): The cobra model object.
        fpkm_df (pd.DataFrame): DataFrame with 'Gene_ID' and 'Expression' columns.

    Returns:
        pd.DataFrame: A DataFrame with reaction IDs as the index and a calculated
                      'g_value' column for use in the objective function.
    """
    print("Starting data pre-processing...")
    start_time = time.time()

    # 1. Efficiently create the FPKM dictionary with the 'G_' prefix.
    # Using a Pandas Series is generally faster for this conversion.
    fpkm_dic = pd.Series(
        fpkm_df["Expression"].values, index='G_' + fpkm_df["Gene_ID"]
    ).to_dict()

    # 2. Cap expression values at the 95th percentile to handle outliers.
    # It's more robust to calculate percentile on non-zero values.
    non_zero_values = [v for v in fpkm_dic.values() if v > 0]
    if non_zero_values:
        cap = np.percentile(non_zero_values, 95)
        print(f"Capping FPKM values at 95th percentile: {cap:.4f}")
        fpkm_dic = {k: min(v, cap) for k, v in fpkm_dic.items()}

    # 3. Calculate 'g' value for each reaction based on its GPR.
    reaction_data = []
    known_g_values = []  # Store g-values for reactions with complete GPR data

    for reaction in model.reactions:
        parsed_rule = parse_gpr(reaction.gene_reaction_rule)

        # If there is no GPR, we cannot calculate 'g' from expression data.
        if not parsed_rule:
            reaction_data.append({'id': reaction.id, 'g_value': np.nan})
            continue

        or_values = []
        is_rule_complete = True
        for and_group in parsed_rule:
            min_val_group = []
            for gene_id in and_group:
                gene_key = f'G_{gene_id}'
                if gene_key in fpkm_dic:
                    min_val_group.append(fpkm_dic[gene_key])
                else:
                    # If any gene in an AND group is missing from FPKM data,
                    # the complex is considered non-functional.
                    is_rule_complete = False
                    break
            
            if not is_rule_complete:
                break
            
            if min_val_group:
                # The value for an AND complex is the minimum expression
                or_values.append(min(min_val_group))

        if is_rule_complete and or_values:
            # The final 'g' value is the sum of the OR'd components (isozymes)
            g_value = sum(or_values)
            reaction_data.append({'id': reaction.id, 'g_value': g_value})
            known_g_values.append(g_value)
        else:
            # Mark as NaN if any gene in the rule was not in the FPKM dictionary
            reaction_data.append({'id': reaction.id, 'g_value': np.nan})

    # 4. Create a DataFrame from the processed data.
    reaction_df = pd.DataFrame(reaction_data).set_index('id')

    # 5. Impute missing 'g' values.
    # Reactions without GPR or with incomplete data get the median 'g' value
    # from all other successfully calculated reactions.
    if known_g_values:
        e_g_median = np.median(known_g_values)
        print(f"Median 'g' value (E_g) for imputation: {e_g_median:.4e}")
        # FIX: Use direct assignment to avoid pandas FutureWarning
        reaction_df['g_value'] = reaction_df['g_value'].fillna(e_g_median)
    else:
        # Fallback if no GPR data could be used at all.
        print("Warning: No GPR data could be mapped. Falling back to g_value=1.0")
        # FIX: Use direct assignment to avoid pandas FutureWarning
        reaction_df['g_value'] = reaction_df['g_value'].fillna(1.0)

    # Add epsilon to all g_values to prevent log(0) in the objective function.
    reaction_df['g_value'] += EPSILON_g

    end_time = time.time()
    print(f"Pre-processing finished in {end_time - start_time:.2f} seconds.")
    return reaction_df


def updateModel(model_default: cobra.Model, mediumFile: str) -> cobra.Model:
    """
    Updates the model by standardizing reaction IDs, bounds, and setting the medium.
    """
    model = model_default.copy()

    # Prefix all reaction IDs with 'R_' if not already present.
    # This standardization helps in matching with other data sources.
    for reaction in model.reactions:
        if not reaction.id.startswith('R_'):
            reaction.id = f'R_{reaction.id}'
        
        # Set bounds to a consistent, wide range for the solver.
        #if reaction.lower_bound < 0:
        #    reaction.lower_bound = REACTION_DEFAULT_LOWER_BOUND
        #if reaction.upper_bound > 0:
        #    reaction.upper_bound = REACTION_DEFAULT_UPPER_BOUND

        if (reaction.lower_bound<0 and reaction.upper_bound>0):
            reaction.bounds = (REACTION_DEFAULT_LOWER_BOUND,REACTION_DEFAULT_UPPER_BOUND)
        if (reaction.lower_bound>=0 and reaction.upper_bound>0):
            reaction.bounds = (0,REACTION_DEFAULT_UPPER_BOUND)
        if (reaction.lower_bound<0 and reaction.upper_bound<=0):
            reaction.bounds = (REACTION_DEFAULT_LOWER_BOUND,0)
    
    # Set culture medium if a medium file is provided.
    if mediumFile != 'NA' and Path(mediumFile).exists():
        print(f"Applying medium from: {mediumFile}")
        medium_df = pd.read_csv(mediumFile, sep="\t", lineterminator='\n')
        
        # The 'medium' property in cobrapy is the recommended way to set exchanges.
        # It automatically handles setting uptake rates. For uptake, the lower bound
        # is set to a negative value.
        #medium_dict = {
        #    f"R_{rxn_id}": abs(REACTION_DEFAULT_LOWER_BOUND)
        #    for rxn_id in medium_df['Reaction_ID']
        #}
        #model.medium = medium_dict
        # Reset bounds of all exchange reactions to (0, 1000)
        for reaction in model.exchanges:
            reaction.bounds = (0, REACTION_DEFAULT_UPPER_BOUND)
        
        # Adjust bounds for reactions specified in the medium file
        for reaction_id in medium_df['Reaction_ID']:
            reaction_id_prefixed = f'R_{reaction_id}'
            if reaction_id_prefixed in model.reactions:
                model.reactions.get_by_id(reaction_id_prefixed).bounds = (REACTION_DEFAULT_LOWER_BOUND, 0)    
    else:
        print("No medium file provided or found. Using model's default exchanges.")
    
    return model


def getPrimalValues(model: cobra.Model,knownFluxes_df: pd.DataFrame) -> List[float]:
    """
    Runs a simple FBA to get an initial guess for the flux values (primal values).
    A good initial guess can significantly speed up the non-linear solver.
    """
    print("Solving with standard FBA to get an initial guess (x0)...")

    modelPV=model.copy()
    objective_reaction_id = list(model.objective.expression.as_coefficients_dict().keys())[0].name
    print(objective_reaction_id)
    # Add constraints for known fluxes, if any.
    if not knownFluxes_df.empty:
        print("Adding known flux constraints...")
        for _, row in knownFluxes_df.iterrows():
            reaction_id = f"R_{row['Reaction_ID']}"
            flux_val = row["Reaction_Flux"]
            print(reaction_id,flux_val)
            if reaction_id == objective_reaction_id:#next(iter(model.objective.variables)).name:
                modelPV.reactions.get_by_id(reaction_id).bounds=(0,flux_val)
            else:
                modelPV.reactions.get_by_id(reaction_id).bounds=(flux_val,flux_val)
    else:
        print(f"Warning: Known flux df is empty.")

    print(modelPV.medium)
    solution = modelPV.optimize()
    if solution.status == "optimal":
        print(f"FBA solution found. Objective value: {solution.objective_value:.4f}. Status: {solution.status}")
        #print(modelPV.summary())
        # Create a list of forward and reverse fluxes based on FBA solution.
        # If flux > 0, vf = flux, vr = 0. If flux < 0, vf = 0, vr = -flux.
        x0 = []
        for n,rxn in enumerate(model.variables):
            x0.append( rxn.primal )          
        return x0
    else:
        print(f"Could not solve FBA for initial guess. Using zero vector. Status: {solution.status}")
        return [0.0] * (len(modelPV.reactions) * 2)



def setVariablesAndObjective(model: cobra.Model, reaction_g_values: pd.DataFrame) -> Tuple:
    """
    Sets up CasADi variables for forward and reverse fluxes and builds the
    entropy-based objective function.
    """
    v_list, v_dic = [], {}
    lbx, ubx = [], []
    g_constraints, lbg, ubg = [], [], []
    
    f_obj = ca.SX(0)  # Initialize objective function as a CasADi symbolic variable.

    for reaction in model.reactions:
        # Retrieve the pre-calculated 'g' value for this reaction.
        ge = reaction_g_values.loc[reaction.id, 'g_value']
        # Set up symbolic variables for forward and reverse fluxes.
        vf = ca.SX.sym(reaction.id)
        vr = ca.SX.sym(reaction.reverse_id)
        
        v_list.extend([vf, vr])
        v_dic[reaction.id] = vf
        v_dic[reaction.reverse_id] = vr
        
        # Set lower and upper bounds for the individual vf and vr variables.
        # These are always non-negative.
        lbx.extend([0.0, 0.0])
        ubx.extend([REACTION_DEFAULT_UPPER_BOUND, REACTION_DEFAULT_UPPER_BOUND])
        
        # Add a constraint for the NET flux (vf - vr) for ALL reactions
        # to enforce the bounds from the cobra model.
        if reaction.lower_bound == 0 and reaction.upper_bound > 0: 
            g_constraints.append(vf - vr)
            lbg.append(0)
            ubg.append(REACTION_DEFAULT_UPPER_BOUND)
        elif reaction.lower_bound < 0 and reaction.upper_bound == 0:
            g_constraints.append(vf - vr)
            lbg.append(REACTION_DEFAULT_LOWER_BOUND)
            ubg.append(0)
        elif reaction.lower_bound == 0 and reaction.upper_bound == 0:
            g_constraints.append(vf - vr)
            lbg.append(0)
            ubg.append(0)            
        
        # Add terms to the objective function: -v * (log(v) - log(g))
        # This is the Shannon entropy term. We use `ca.if_else` to avoid log(0).
        f_obj += ca.if_else(vf > EPSILON, -vf * (ca.log(vf) - ca.log(ge)), EPSILON)
        f_obj += ca.if_else(vr > EPSILON, -vr * (ca.log(vr) - ca.log(ge)), EPSILON)

    return v_list, v_dic, lbx, ubx, f_obj, g_constraints, lbg, ubg


def createConstraints(
    model: cobra.Model,
    v_dic: Dict[str, ca.SX],
    g_constraints: List[ca.SX],
    lbg: List[float],
    ubg: List[float],
    knownFluxes_df: pd.DataFrame
) -> Tuple[List[float], List[float], List[ca.SX]]:
    """
    Creates stoichiometric and known flux constraints for the NLP problem.
    This version builds constraints directly from model properties for efficiency.
    """
    print("Creating stoichiometric constraints...")
    # Stoichiometric constraints (S*v = 0)
    for met in model.metabolites:
        # FIX: Correctly iterate over the reactions a metabolite participates in.
        # `met.reactions` is a frozenset of Reaction objects.
        participating_reactions = met.reactions
        
        constraint_expr = ca.SX(0)
        for rxn in participating_reactions:
            # Get the stoichiometric coefficient of the current metabolite in this reaction.
            coeff = rxn.metabolites[met]
            # Net flux is the difference between forward and reverse variables.
            net_flux = v_dic[rxn.id] - v_dic[rxn.reverse_id]
            constraint_expr += coeff * net_flux
            
        g_constraints.append(constraint_expr)
        lbg.append(0.0)  # Steady-state assumption
        ubg.append(0.0)

    # Add constraints for known fluxes, if any.
    if not knownFluxes_df.empty:
        print("Adding known flux constraints...")
        for _, row in knownFluxes_df.iterrows():
            reaction_id = f"R_{row['Reaction_ID']}"
            flux_val = row["Reaction_Flux"]
            
            if reaction_id in v_dic:
                v_f = v_dic[model.reactions.get_by_id(reaction_id).id]
                v_r = v_dic[model.reactions.get_by_id(reaction_id).reverse_id]
                new_constraint = (v_f - v_r) - flux_val
                g_constraints.append(new_constraint)
                lbg.append(0.0) # The difference should be exactly 0.
                ubg.append(0.0)
            else:
                print(f"Warning: Known flux Reaction_ID '{row['Reaction_ID']}' not in model.")

    return ubg, lbg, g_constraints

def optPheFlux(
    model: cobra.Model,
    reaction_g_values: pd.DataFrame,
    init_time: float,
    knownFluxes_df: pd.DataFrame,
    ipoptParams: dict
) -> Tuple:
    """
    Main optimization function that sets up and solves the non-linear problem.
    """
    # 1. Set up variables and objective function
    v, v_dic, lbx, ubx, f_obj, g, lbg, ubg = setVariablesAndObjective(model, reaction_g_values)
    ncbmb=len(lbg) # number of constraints before mass balances
    
    # 2. Create stoichiometric and other constraints
    ubg, lbg, g = createConstraints(model, v_dic, g, lbg, ubg, knownFluxes_df)

    # 3. Define the Non-Linear Programming (NLP) problem for CasADi
    nlp = {
        'x': ca.vertcat(*v),    # Decision variables
        'f': -f_obj,            # Objective function (we maximize, so we negate)
        'g': ca.vertcat(*g),    # Constraints
    }
    # 4. Configure the IPOPT solver
    opts = {
        'ipopt': ipoptParams,
        'print_time': False
    }
    solver = ca.nlpsol('solver', 'ipopt', nlp, opts)
    
    # 5. Get initial guess (x0) from a standard FBA
    x0 = getPrimalValues(model,knownFluxes_df)

    # 6. Solve the problem
    print("## Solving the non-linear optimization problem with IPOPT...")
    start_opt_time = time.time()
    sol = solver(x0=x0, ubg=ubg, lbg=lbg, lbx=lbx, ubx=ubx)
    end_opt_time = time.time()
    
    optimization_time = end_opt_time - start_opt_time
    total_time = end_opt_time - init_time
    status = solver.stats()['return_status']
    success = solver.stats()['success']
    
    print(f"Solver finished. Status: {status} | Success: {success}")
    print(f"Optimization took: {optimization_time:.2f} s")
    
    # 7. Process and return results
    predicted_fluxes = sol['x']
    
    forward_reverse_fluxes = pd.DataFrame(columns=["Reaction_ID","forward", "reverse"])
    flux_results = {}
    for i in range(0, len(v), 2):
        reaction_name = str(v[i])
        # Net flux = forward flux - reverse flux
        net_flux = float(predicted_fluxes[i] - predicted_fluxes[i+1])
        flux_results[reaction_name] = net_flux
        forward_reverse_fluxes.loc[reaction_name] = [reaction_name,float(predicted_fluxes[i]), float(predicted_fluxes[i+1])]
        
    flux_series = pd.Series(flux_results)

    lam_g = sol['lam_g'] # Lagrange multipliers of constraints
    lagrange= pd.DataFrame(columns=["Lagrange_coef"])
    for i,metabolite in enumerate(model.metabolites):
        metabolite_name = metabolite.id 
        lagrange_value  = float(lam_g[i+ncbmb]) 
        lagrange.loc[metabolite_name] = [lagrange_value]

    return flux_series, optimization_time, total_time, status, success, lbx, ubx, forward_reverse_fluxes,lagrange

def recordTable(record: pd.DataFrame, condition: str, lbx: list, ubx: list, time: float, status: str, memory_mb: float) -> pd.DataFrame: 
    """Records summary statistics for a single run."""
    # Count variables that are not fixed (lower bound != upper bound)
    num_variables = sum(1 for i in range(len(lbx)) if lbx[i] != ubx[i])
    
    new_row = pd.DataFrame({
        'Condition': [condition],
        'N° variables': [num_variables],
        'Time (s)': [time],
        'Status': [status],
        'Memory (MB)': [memory_mb]  # Add the new column here
    })
    
    return pd.concat([record, new_row], ignore_index=True)


def actuallyTime() -> str:
    """Returns the current timestamp in a formatted string for logging."""
    return datetime.now().strftime("[%Y/%m/%d %H:%M:%S]")


def getFluxes(inputFileName: str, processDir: str, prefix: str, ipoptParams: dict=None) -> None:
    """
    Main orchestration function that reads the input file and processes each condition.
    """
    if ipoptParams is None: 
        ipoptParams =   {'print_level': 0, 'sb': 'yes'}       
    processStart = time.time()
    summary_record = pd.DataFrame()
    
    # Load the main input file describing all conditions to run.
    try:
        inputData = pd.read_csv(inputFileName, sep="\t", lineterminator='\n', na_filter=False)
    except FileNotFoundError:
        print(f"Error: Input file not found at {inputFileName}")
        return

    opt_times, total_times = [], []
    allFluxes,allForRev,allLagran = [],[],[] 

    print(f"Found {len(inputData)} conditions to process.")
    for i, row in inputData.iterrows():
        # --- Load information for the current condition ---
        condition = row["Condition"]
        geneExpFile = row["GeneExpFile"]
        mediumFile = row["Medium"]
        networkFile = row["Network"]
        organism = row["Organism"]
        knownFluxesFile = row["KnownFluxes"]

        print('\n' + f' STARTING: {organism} - {condition} '.center(80, '='))
        
        # --- Load Model and Data ---
        print(f"{actuallyTime()} Loading metabolic model: {Path(networkFile).name}")
        model_default = cobra.io.read_sbml_model(networkFile)
        
        print(f"{actuallyTime()} Loading transcriptomic data: {Path(geneExpFile).name}")
        fpkm_df = pd.read_csv(geneExpFile, sep="\t", lineterminator='\n')
        
        knownFluxes_df = pd.DataFrame()
        if knownFluxesFile != 'NA' and Path(knownFluxesFile).exists():
            print(f"{actuallyTime()} Loading known fluxes: {Path(knownFluxesFile).name}")
            knownFluxes_df = pd.read_csv(knownFluxesFile, sep="\t", lineterminator='\n')
        
        init_time = time.time()

        # --- Update and Pre-process ---
        print(f"{actuallyTime()} Updating model (bounds, medium)...")
        model = updateModel(model_default, mediumFile)
        
        # THIS IS THE KEY OPTIMIZATION:
        reaction_g_values = preprocess_reaction_data(model, fpkm_df)
        del fpkm_df, model_default # Free up memory
        gc.collect()
        
        # --- Run Flux Prediction ---
        print(f"{actuallyTime()} Running TeraFlux optimization...")
        fluxes, opt_time, t_time, status, success, lbx, ubx,forward_reverse_fluxes,lagrange = optPheFlux(
            model, reaction_g_values, init_time, knownFluxes_df,ipoptParams
        )

        # --- Save Results ---
        print(f"{actuallyTime()} Saving results...")
        output_dir = Path(processDir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save fluxes
        allFluxes.append(fluxes)
        allForRev.append(forward_reverse_fluxes)
        allLagran.append(lagrange)
        
        # Update summary table
        # Get memory usage at the end of the iteration
        process = psutil.Process(os.getpid())
        memory_mb = process.memory_info().rss / (1024 * 1024)

        # Update summary table with memory usage
        summary_record = recordTable(summary_record, f"{organism}_{condition}", lbx, ubx, t_time, status, memory_mb)       
        
        opt_times.append(opt_time)
        total_times.append(t_time)
        print(f"{actuallyTime()} {organism} - {condition} processed.")
        print('=' * 80)

    # --- Final Summary ---
    print('\n' + ' PROCESS COMPLETE '.center(80, '='))
    if not summary_record.empty:
        recordFile = Path(processDir) / f"{prefix}.log.tsv"
        summary_record.to_csv(recordFile, sep='\t', index=False)
        print(f"Summary log saved to: {recordFile}")

    processFinish = time.time()
    processTime = processFinish - processStart
    
    if total_times:
        print(f'Average time per optimization: {np.mean(opt_times):.2f} s')
        print(f'Average time per condition: {np.mean(total_times):.2f} s')
    print(f'Total process time: {processTime / 60:.2f} min ({processTime / 3600:.2f} h)')
    
    mem_usage = psutil.Process(os.getpid()).memory_info().rss / (1024 * 1024)
    print(f"Final memory usage: {mem_usage:.2f} MB")

    # Concatenate all DataFrames in the list into a single DataFrame
    allFluxes = pd.concat(allFluxes, axis=1)  
    allForRev = pd.concat(allForRev, axis=1)  
    allLagran = pd.concat(allLagran, axis=1)  
    # Transpose the DataFrame to have reactions as rows and samples as columns
    #allFluxes = allFluxes.T 
    # Save results   
    allFluxes.to_csv(output_dir / f'{prefix}.fluxes.csv' ) 
    allForRev.to_csv(output_dir / f'{prefix}.forward_reverse_fluxes.csv' ) 
    allLagran.to_csv(output_dir / f'{prefix}.Lagrange.csv' ) 

    
    return(allFluxes,allForRev,allLagran)

def main():
    """Main function to run the script from the command line."""
    parser = argparse.ArgumentParser(
        description="Run TeraFlux to predict metabolic fluxes from transcriptomic data."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the main input file listing all conditions to process."
    )
    parser.add_argument(
        "-o", "--output_dir",
        type=str,
        default="teraflux_results",
        help="Directory to save the output files."
    )
    parser.add_argument(
        "-p", "--prefix",
        type=str,
        default="teraflux_run",
        help="Prefix for the summary log file."
    )
    parser.add_argument(
        "-i", "--ipoptParams",
        type=dict,
        help="Dictionary of IPOPT parameters."
    )
    args = parser.parse_args()

    print('Welcome to Optimized TeraFlux!')
    getFluxes(args.input_file, args.output_dir, args.prefix, args.ipoptParams)


if __name__ == "__main__":
    main()
