# Global Imports
import argparse
import os
import sys
import yaml
import random
import logging
import hashlib


def main(args):
    input_contrast_tsv = args.input_contrast_tsv
    outdir_path = args.outdir_path
    name = args.name
    skip_tf_activity = args.skip_tf_activity
    skip_pathway_activity = args.skip_pathway_activity
    skip_gsea = args.skip_gsea
    skip_ora = args.skip_ora
    tf_annotation_resource = "collectri"
    tf_method = "ulm"
    tf_additional_args = {}
    pathway_annotation_resource = "progeny"
    pathway_method = "mlm"
    pathway_additional_args = {"top": 500}
    gsea_resource = "msigdb"
    gsea_msigdb_sets = [
        "hallmark",
        "kegg_pathways", 
        "go_biological_process", 
        "go_molecular_function", 
        "go_cellular_component",
        "reactome_pathways",
    ]
    gsea_additional_args = {"padj": 0.05, "top_n": None}
    ora_resource = "msigdb"
    ora_msigdb_sets = [
        "hallmark",
        "kegg_pathways", 
        "go_biological_process", 
        "go_molecular_function", 
        "go_cellular_component",
        "reactome_pathways",
    ]
    ora_additional_args = {"padj": 0.05, "top_n": None}

    # Make the output directory if it doesn't exist
    if not os.path.exists(outdir_path):
        os.makedirs(outdir_path)

    # Log file
    if os.path.exists(os.path.join(outdir_path, "functional_analysis.log")):
        os.remove(os.path.join(outdir_path, "functional_analysis.log"))

    # Make a workflow hash
    workflow_hash = hashlib.sha1()
    workflow_hash.update(str(random.getrandbits(128)).encode("utf-8"))
    
    # Log file
    if os.path.exists(os.path.join(outdir_path, "functional_analysis.log")):
        os.remove(os.path.join(outdir_path, "functional_analysis.log"))
    logging.basicConfig(
        filename=os.path.join(outdir_path, "functional_analysis.log"),
        level=logging.INFO, 
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    logging.info(f"Workflow hash: {workflow_hash.hexdigest()}")

    # Make and log params
    import scanpy as sc
    import decoupler as dc
    data_params = {
        "input_contrast_tsv": input_contrast_tsv,
        "outdir_path": outdir_path,
    }
    tf_activity_params = {
        "tf_annotation_resource": tf_annotation_resource,
        "tf_method": tf_method,
        "tf_additional_args": tf_additional_args,
    }
    pathway_activity_params = {
        "pathway_annotation_resource": pathway_annotation_resource,
        "pathway_method": pathway_method,
        "pathway_additional_args": pathway_additional_args,
    }
    gsea_params = {
        "gsea_resource": gsea_resource,
        "gsea_msigdb_sets": gsea_msigdb_sets,
        "gsea_additional_args": gsea_additional_args,
    }
    ora_params = {
        "ora_resource": ora_resource,
        "ora_msigdb_set": ora_msigdb_sets,
        "ora_additional_args": ora_additional_args,
    }
    version_params = {
        "worklow_hash": workflow_hash.hexdigest(),
        "Python": sys.version[:5],
        "Scanpy": sc.__version__,
        "decoupler": dc.__version__,
    }
    params = {
        "data": data_params, 
        "tf_activity": tf_activity_params if not skip_tf_activity else None,
        "pathway_activity": pathway_activity_params if not skip_pathway_activity else None,
        "gsea": gsea_params if not skip_gsea else None,
        "ora": ora_params if not skip_ora else None,
        "version": version_params
    }
    if not os.path.exists(os.path.join(outdir_path, "functional_analysis.yaml")):
        logging.info("Writing params to {}".format(outdir_path))
        with open(
            os.path.join(outdir_path, "functional_analysis.yaml"), "w"
        ) as outfile:
            yaml.dump(params, outfile, default_flow_style=False)
    else:
        logging.info("params.yaml already exists, will not overwrite")

    # The data to load in formatted as a 10x directory, change if you have a different format for reading in data
    import pandas as pd
    logging.info(f"Reading in {input_contrast_tsv} as a dataframe")
    results_df = pd.read_csv(input_contrast_tsv, sep="\t", index_col=0)
    results_df.head()
    mat = results_df[['stat']].T.rename(index={'stat': name})

    if not skip_tf_activity:
        logging.info("Running TF activity analysis")
        
        # Retrieve GRN
        logging.info(f"Retrieving {tf_annotation_resource} GRN")
        if tf_annotation_resource == "dorothea":
            net = dc.get_dorothea(organism='human', split_complexes=False)
        elif tf_annotation_resource == "collectri":
            net = dc.get_collectri(organism='human', split_complexes=False)
        else:
            raise ValueError(f"tf_annotation_resource {tf_annotation_resource} not supported")
        
        # Infer pathway activities with ulm
        logging.info(f"Running {tf_method} on {tf_annotation_resource} GRN")
        tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=net)
        
        # Plot the TF activity barplot
        logging.info(f"Plotting TF activity barplot")
        dc.plot_barplot(tf_acts, name, top=25, vertical=True, save=os.path.join(outdir_path, f"{name}_tf_activity_barplot.png"))

        # Save the TF activity table
        tf_acts.T.to_csv(os.path.join(outdir_path, f"{name}_tf_activity.tsv"), sep="\t")

    if not skip_pathway_activity:
        logging.info("Running pathway activity analysis")

        # Retrieve weights
        if pathway_annotation_resource == "progeny":
            logging.info(f"Retrieving {pathway_annotation_resource} weights")
            net = dc.get_progeny(top=pathway_additional_args["top"])
        else:
            raise ValueError(f"pathway_annotation_resource {pathway_annotation_resource} not supported")

        # Infer pathway activities with mlm
        logging.info(f"Running {pathway_method} on {pathway_annotation_resource} weights")
        pathway_acts, pathway_pvals = dc.run_mlm(mat=mat, net=net)

        # Plot a barplot of the pathway activities
        logging.info(f"Plotting pathway activity barplot")
        dc.plot_barplot(pathway_acts, name, top=25, vertical=False, save=os.path.join(outdir_path, f"{name}_pathway_activity_barplot.png"))

        # Save the pathway activities to a file
        pathway_acts.T.to_csv(os.path.join(outdir_path, f"{name}_pathway_activity.tsv"), sep="\t")

    if not skip_gsea:
        logging.info("Running GSEA")
        for gsea_msigdb_set in gsea_msigdb_sets:
            logging.info(f"Running GSEA on {gsea_msigdb_set}")
            # Retrieve resource
            if gsea_resource == "msigdb":
                logging.info(f"Retrieving {gsea_resource} gene sets")
                gene_sets = dc.get_resource('MSigDB')
                
                # Filter by kegg_pathways
                logging.info(f"Filtering {gsea_resource} gene sets by {gsea_msigdb_set}")
                gene_sets = gene_sets[gene_sets['collection'] == gsea_msigdb_set]

                # Remove duplicated entries
                gene_sets = gene_sets[~gene_sets.duplicated(['geneset', 'genesymbol'])]

                # Rename
                if gsea_msigdb_set == "hallmark":
                    gene_sets.loc[:, 'geneset'] = [name.split('HALLMARK_')[1] for name in gene_sets['geneset']]

            else:
                raise ValueError(f"gsea_resource {gsea_resource} not supported")

            # Choose genes to run GSEA on
            if gsea_additional_args["padj"] is None:
                assert gsea_additional_args["top_n"] is not None
                top_genes = results_df.sort_values("padj").iloc[:gsea_additional_args["top_n"]]
                logging.info(f"Running GSEA on top {gsea_additional_args['top_n']} genes")
            else:
                top_genes = results_df[results_df['padj'] < 0.05]
                logging.info(f"Running GSEA on all genes with padj < {gsea_additional_args['padj']}")
            logging.info(f"Number of genes to run GSEA on: {top_genes.shape[0]}")

            # Run GSEA
            enr_pvals = dc.get_gsea_df(
                df=top_genes,
                stat="log2FoldChange",
                net=gene_sets,
                source='geneset',
                target='genesymbol',
            )

            # Save the enrichment table
            enr_pvals.to_csv(os.path.join(outdir_path, f"{name}_{gsea_resource}_{gsea_msigdb_set}_gsea_enrichment.tsv"), sep="\t")
            
            # Keep top 25 based on NES
            enr_pvals = enr_pvals.sort_values("NES", ascending=False).iloc[:25]
            enr_pvals["size"] = 1

            # Plot the dotplot
            dc.plot_dotplot(
                enr_pvals, 
                x='NES', 
                y='Term',
                s="size",
                c='FDR p-value', 
                scale = 0.5, 
                figsize=(7,14), 
                save=os.path.join(outdir_path, f"{name}_{gsea_resource}_{gsea_msigdb_set}_gsea_dotplot.png")
            )

    if not skip_ora:
        logging.info("Running ORA")

        for ora_msigdb_set in ora_msigdb_sets:

            logging.info(f"Running ORA on {ora_msigdb_set}")
            
            # Retrieve resource
            if ora_resource == "msigdb":

                # Retrieve resource
                logging.info(f"Retrieving {ora_resource} gene sets")
                gene_sets = dc.get_resource('MSigDB')
                
                # Filter by kegg_pathways
                logging.info(f"Filtering {ora_resource} gene sets by {ora_msigdb_set}")
                gene_sets = gene_sets[gene_sets['collection'] == ora_msigdb_set]

                # Remove duplicated entries
                gene_sets = gene_sets[~gene_sets.duplicated(['geneset', 'genesymbol'])]

                # Rename
                if ora_msigdb_set == "hallmark":
                    gene_sets.loc[:, 'geneset'] = [name.split('HALLMARK_')[1] for name in gene_sets['geneset']]

            else:
                raise ValueError(f"ora_resource {ora_resource} not supported")
            
            # Choose genes to run ORA on
            if ora_additional_args["padj"] is None:
                assert ora_additional_args["top_n"] is not None
                top_genes = results_df.sort_values("padj").iloc[:ora_additional_args["top_n"]]
                logging.info(f"Running ORA on top {ora_additional_args['top_n']} genes")
            else:
                top_genes = results_df[results_df['padj'] < 0.05]
                logging.info(f"Running ORA on all genes with padj < {ora_additional_args['padj']}")
            logging.info(f"Number of genes to run ORA on: {top_genes.shape[0]}")
            
            # Run ora
            enr_pvals = dc.get_ora_df(
                df=top_genes,
                net=gene_sets,
                source='geneset',
                target='genesymbol'
            )

            # Save the enrichment table
            enr_pvals.to_csv(os.path.join(outdir_path, f"{name}_{ora_resource}_{ora_msigdb_set}_ora_enrichment.tsv"), sep="\t")

            # Keep top 25 based on combined score
            enr_pvals = enr_pvals.sort_values("Combined score", ascending=False).iloc[:25]

            # Plot the dotplot
            dc.plot_dotplot(
                enr_pvals, 
                x='Combined score', 
                y = 'Term', 
                s='Odds ratio', 
                c = 'FDR p-value', 
                scale = 0.5, 
                figsize=(7,14), 
                save=os.path.join(outdir_path, f"{name}_{ora_resource}_{ora_msigdb_set}_ora_dotplot.png")
            )

    logging.info("Done with functional analysis")


if __name__ == "__main__":
    # Setting up argparse
    parser = argparse.ArgumentParser(
        description="Running functional analysis on provided contrast."
    )
    parser.add_argument(
        "--input_contrast_tsv", required=True, help="Path to the input contrast tsv file."
    )
    parser.add_argument(
        "--outdir_path", required=True, help="Path to the output directory."
    )
    parser.add_argument(
        "--name", required=True, help="Name of the contrast."
    )
    parser.add_argument(
        "--skip_tf_activity", required=False, default=False, help="Boolean specifying whether to skip TF activity analysis."
    )
    parser.add_argument(
        "--skip_pathway_activity", required=False, default=False, help="Boolean specifying whether to skip pathway activity analysis."
    )
    parser.add_argument(
        "--skip_gsea", required=False, default=False, help="Boolean specifying whether to skip GSEA."
    )
    parser.add_argument(
        "--skip_ora", required=False, default=False, help="Boolean specifying whether to skip ORA."
    )
    args = parser.parse_args()
    main(args)
    