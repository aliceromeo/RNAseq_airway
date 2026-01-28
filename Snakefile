import pandas as pd
# Load parameters
configfile: "config/config.yaml"
# Get samples IDs
metadata = pd.read_csv(config.get("metadata", "metadata.csv"))
SAMPLES = metadata['Run (SRA ID)'].tolist()

# Counts matrices for each sample
rule all:
    input:
        # Salmon quantification
        expand("data/salmon_results/{sample}_quant/quant.sf", sample=SAMPLES),
        
        # Quality control (MultiQC)
        "reports/multiqc_report.html",
        
        # Data import and DESeq
        "data/txi_output.rds",
        "data/dds_results.rds",
        "results/tables/complete_results_DGE.csv",
        expand("results/tables/significant_genes_FC{fc}.csv", fc=config["fc_thresholds"]),
        
        # EDA
        "results/plots/2.plotPCA.png",
        "results/plots/1.screePlot.png",
        "results/plots/3.pcaPairs.png",
        "results/plots/4.combinedPlotLoad.png",
        
        # Functional enrichment (GO e GSEA)
        "results/plots/1.go_dotplot.png",
        "results/plots/4.gsea_top_pathways.png",
        "results/plots/2.gsea_dotplot.png",
        "results/plots/3.gsea_ridgeplot.png",
        "results/plots/5.heatmap_core.png"

        
# RULE 1: download fastq files
rule download_fastq:
    output:
        r1="data/fastq_raw/{sample}_1.fastq.gz",
        r2="data/fastq_raw/{sample}_2.fastq.gz"
    threads: config["threads_default"]
    conda: "envs/preproc_env.yaml"
    shell:
        """
        parallel-fastq-dump --sra-id {wildcards.sample} \
            -t {threads} --outdir data/fastq_raw/ --split-files --gzip
        """

# RULE 2: Download references Genome e Transcriptome
rule download_refs:
    output:
        genome="data/references/GRCh38.primary_assembly.genome.fa.gz",
        transcripts="data/references/gencode.v44.transcripts.fa.gz",
        gtf_gz="data/references/gencode.v44.annotation.gtf.gz",
        gtf="data/references/gencode.v44.annotation.gtf"
    conda: "envs/preproc_env.yaml" 
    shell:
        """
        wget -c -P data/references/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh38.primary_assembly.genome.fa.gz
        wget -c -P data/references/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz
        wget -c -P data/references/ https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz
        
        # Unzip GTF file
        gunzip -c {output.gtf_gz} > {output.gtf}
        """

# RULE 3: Creating Salmon index
rule salmon_index:
    input:
        genome="data/references/GRCh38.primary_assembly.genome.fa.gz",
        transcripts="data/references/gencode.v44.transcripts.fa.gz"
    output:
        directory(config["index_dir"])
    threads: config["threads_default"]
    conda: "envs/preproc_env.yaml"
    shell:
        """
        gunzip -c {input.genome} > data/references/genome.fa
        gunzip -c {input.transcripts} > data/references/transcripts.fa
        
        grep ">" data/references/genome.fa | cut -d " " -f 1 | sed 's/>//g' > data/decoys.txt
        cat data/references/transcripts.fa data/references/genome.fa > data/references/transcriptome_genome.fa
        
        salmon index -t data/references/transcriptome_genome.fa -d data/decoys.txt -i {output} -p {threads} --gencode
        
        # Remove unzipped and temporary files
        rm data/references/genome.fa data/references/transcripts.fa data/references/transcriptome_genome.fa
        """

# RULE 4: Quantification
rule salmon_quant:
    input:
        index=config["index_dir"],
        r1="data/fastq_raw/{sample}_1.fastq.gz",
        r2="data/fastq_raw/{sample}_2.fastq.gz"
    conda: "envs/preproc_env.yaml" 
    output:
        sf="data/salmon_results/{sample}_quant/quant.sf",
        meta="data/salmon_results/{sample}_quant/lib_format_counts.json"
    threads: config["threads_default"]
    shell:
        """
        salmon quant -i {input.index} -l A \
            -1 {input.r1} -2 {input.r2} \
            --validateMappings --gcBias --seqBias \
            --threads {threads} -o data/salmon_results/{wildcards.sample}_quant
        """

# RULE 5: FASTQC 
rule fastqc:
    input:
        r1="data/fastq_raw/{sample}_1.fastq.gz",
        r2="data/fastq_raw/{sample}_2.fastq.gz"
    output:
        html1="reports/fastqc/{sample}_1_fastqc.html",
        zip1="reports/fastqc/{sample}_1_fastqc.zip",
        html2="reports/fastqc/{sample}_2_fastqc.html",
        zip2="reports/fastqc/{sample}_2_fastqc.zip"
    conda: "envs/preproc_env.yaml"
    shell:
        "fastqc {input.r1} {input.r2} -o reports/fastqc/"
        
# RULE 6: Report aggregation (MultiQC)
rule multiqc:
    input:
        # All fastqc and salmon reports
        expand("reports/fastqc/{sample}_{n}_fastqc.zip", sample=SAMPLES, n=[1,2]),
        expand("data/salmon_results/{sample}_quant/quant.sf", sample=SAMPLES)
    output:
        report("reports/multiqc_report.html", caption="captions/multiqc.rst", category="1. Quality Control")
    conda: "envs/preproc_env.yaml"
    shell:
        "multiqc reports/fastqc/ data/salmon_results/ -o reports/ -n multiqc_report.html -f"
        
# Starting analysis in R  
# RULE 7: Data import in R
rule r_data_import:
    input:
        salmon_files = expand("data/salmon_results/{sample}_quant/quant.sf", sample=SAMPLES),
        samples_meta = "data/samples.csv",
        gtf = "data/references/gencode.v44.annotation.gtf"
    output:
        rds = "data/txi_output.rds"
    conda: "envs/r_env.yaml"
    script: "scripts/01_data_import.R"

# RULE 8: DESeq2
rule r_deseq:
    input:
        rds_import = "data/txi_output.rds"
    output:
        dds_final = "data/dds_results.rds",
        csv_all = "results/tables/complete_results_DGE.csv",
        csv_sign = expand("results/tables/significant_genes_FC{fc}.csv", fc=config["fc_thresholds"]),
        # add report
        pca_plot = report("results/plots/2.plotPCA.png",
        		  caption="captions/pca.rst", 
        		  category = "2. Differential Expression"),
        scree_plot = report("results/plots/1.screePlot.png",
        		     caption = "captions/scree.rst",
        		     category = "2. Differential Expression"),
        pca_pairs = report("results/plots/3.pcaPairs.png",
        		   caption = "captions/pcapairs.rst",
        		   category = "2. Differential Expression"),
        combined_plot_load = report("results/plots/4.combinedPlotLoad.png",
        		   caption = "captions/loadings.rst",
        		   category = "2. Differential Expression")
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/02_differential_expression.R"

# RULE 9: Functional enrichment (GO e GSEA)
rule r_functional_enrich:
    input:
        dds_results = "data/dds_results.rds"
    output:
        # --- Output GO ---
        go_plot = report("results/plots/1.go_dotplot.png",
        		 caption = "captions/goDot.rst",
         		 category = "3. Functional enrichment"),
         		 
        # --- Output GSEA ---
        gsea_plot = report("results/plots/4.gsea_top_pathways.png",
        		   caption = "captions/gseaFirstPath.rst",
        		   category = "3. Functional enrichment"),
        gsea_dotplot = report("results/plots/2.gsea_dotplot.png",
        		      caption = "captions/gseaDot.rst",
        		      category = "3. Functional enrichment"),
        gsea_ridge = report("results/plots/3.gsea_ridgeplot.png",
        		    caption = "captions/gseaRidge.rst",
        		    category = "3. Functional enrichment"),
        heatmap_core = report("results/plots/5.heatmap_core.png",
        		      caption = "captions/heatmap.rst",
        		      category = "3. Functional enrichment")
    conda:
        "envs/r_env.yaml"
    script:
        "scripts/03_functional_enr.R"
