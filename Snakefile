import os

configfile: "config.yaml"

if os.path.exists("Snakefile_versioned.sk"):
    include: "Snakefile_versioned.sk"

rule sort:
    input:
        "best_param_3plex/{lncRNA}.neg_pos_rand.bed"
    output:
        "best_param_3plex/{lncRNA}.neg_pos_rand.sorted.bed"
    shell: """
        sort -k1,1 -k2,2n {input} | bedtools merge -c 4,5,6 -o first,distinct,distinct > {output}
    """
        
rule shape_mean:
    input:
        shape="best_param_3plex/{species}.{shape}.bw",
        bed="best_param_3plex/{lncRNA}.neg_pos_rand.sorted.bed"
    output:
        "best_param_3plex/{lncRNA}.{species}.{shape}.wig.mean"
    shell: """
        bedtools intersect -sorted -a {input.bed} -b <(bigWigToBedGraph {input.shape} /dev/stdout) -loj \
        | tr -d '\r' \
        | bawk '{{print $1,$2,$3,$4,$6,$10}}' \
        | bedtools merge -c 4,5,6 -o distinct,distinct,mean \
        | sort -k5,5r -k1,1 -k2,2n > {output}
    """
    
rule wig2bw:
    input:
        shape="{prefix}mm10.{suffix}.wig",
        chrom_size="mm10.chrom.sizes"
    output:
        "{prefix}mm10.{suffix}.bw"
    shell: """
        wigToBigWig {input.shape} {input.chrom_size} {output}
    """
    
rule get_all_bw:
    input:
      "best_param_3plex/mm10.{shape}.bw", shape=config["SHAPEFEATURES"]
       
rule aggregate:
    input:
        expand("best_param_3plex/{h_lncRNA}.hg38.{shape}.wig.mean", shape=config["SHAPEFEATURES"], h_lncRNA=config["HUMAN_LNCRNAs"]),
        expand("best_param_3plex/{m_lncRNA}.mm10.{shape}.wig.mean", shape=config["SHAPEFEATURES"], m_lncRNA=config["MOUSE_LNCRNAs"])
    output:
        "best_param_3plex/ALL_shape.aggregated.gz"
    shell:"""
        matrix_reduce -t '*.wig.mean' -l '{input}' | sed 's/.mm10./:/' | sed 's/.hg38./:/' | tr "/" "\t" | tr ":" "\t" \
        | bawk 'BEGIN {{print "lncRNA","Shape","ID","pos_neg","Mean"}} {{print $2,$3,$7,$8,$9}}' | gzip > {output}
    """

rule shape_prep4join:
    input:
        "best_param_3plex/ALL_shape.aggregated.gz"
    output:
        "best_param_3plex/ALL_shape.aggregated.rearranged.gz"   
    shell:"""
        zcat {input} | grep -v , | bawk '{{print $3,$1,$2,$4,$5}}' | sort -k1 | gzip > {output}
    """    
 
#Inputs have different formats, so two rules are needed to correctly process them

#for lncRNA1: zcat best_param_3plex/results_3plex/{lncRNA}_ssmasked-{lncRNA}.neg_pos_rand.bed.tpx.summary.add_zeros.gz \
#             | tr ":" "\t" | bawk '{{print $1,$15,$18}}' > {lncRNA}_stability.summary
#lncRNA1=["7SK", "AC018781.1", "AC087482.1", "AL109615.3", "Bloodlinc", "CPEB2-DT", "Eprn", "HAND2-AS1", "MALAT1", "MIR503HG", "Rn7sk", "SRA1", "Tug1"]

#for lncRNA2: zcat best_param_3plex/results_3plex/{lncRNA}_ssmasked-{lncRNA}.neg_pos_rand.bed.tpx.summary.add_zeros.gz \
#             | bawk '{{print $1,$15,$18}}' > {lncRNA}_stability.summary 
#lncRNA2=["CDKN2B-AS1", "HOTAIR", "lncSmad7", "MEG3", "Meg3", "NEAT1", "TERC"]  

rule stability_prep4join:
    input:
        #expand("best_param_3plex/results_3plex/{h_lncRNA}_ssmasked-{h_lncRNA}.neg_pos_rand.bed.tpx.summary.add_zeros", h_lncRNA=config["HUMAN_LNCRNAs"]),
        #expand("best_param_3plex/results_3plex/{m_lncRNA}_ssmasked-{m_lncRNA}.neg_pos_rand.bed.tpx.summary.add_zeros", m_lncRNA=config["MOUSE_LNCRNAs"])
        expand("best_param_3plex/results_3plex/{h_lncRNA}_stability.summary", h_lncRNA=config["HUMAN_LNCRNAs"]),
        expand("best_param_3plex/results_3plex/{m_lncRNA}_stability.summary", m_lncRNA=config["MOUSE_LNCRNAs"])
    output:
        "best_param_3plex/ALL_stability.aggregated.rearranged.gz" 
    shell:"""
        matrix_reduce '{input}' | grep -v ">" | sort -k1 | bawk 'BEGIN {{print "ID","Stability_best","Stability_norm"}} {{print $1,$2,$3}}' | gzip > {output}
    """

rule join_files:
    input:
        shapes="best_param_3plex/ALL_shape.aggregated.rearranged.gz",
        stability="best_param_3plex/ALL_stability.aggregated.rearranged.gz" 
    output:
        "best_param_3plex/ALL_shape.3plex_stability.gz"
    shell:"""
        join -j 1 <(zcat {input.shapes} | sort) <(zcat {input.stability}| sort) |tr " " "\t" | sort -k2,2 -k3,3 -k4,4r \
        | bawk 'BEGIN {{print "lncRNA","Shape","pos_neg","Mean", "Stability_best","Stability_norm", "ID"}} {{print $2,$3,$4,$5,$6,$7,$1}}'| gzip > {output}
    """

rule tab2matrix:
    input:
        "best_param_3plex/ALL_shape.3plex_stability.gz"
    output:
        "best_param_3plex/ALL_shape.3plex_stability.matrix.gz"
    shell:"""
        zcat {input} | bawk '$1!="lncRNA" && $1!="Shape"  {{print $7","$6","$5","$3","$1,$2,$4}}' \
        | tab2matrix -r "ID,Stability_norm,Stability_best,pos_neg,lncRNA" | tr "," "\t" \
        | bawk '{{print $1,$5,$4,$2,$3,$6~9}}'\
        | gzip > {output}
    """


#Duplicate rule to retrieve the correct genome
rule h_geneName_to_fasta:
        input:
            path="/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/transcripts.longest.fa.gz"
        output:
            "best_param_3plex/hg38.{h_lncRNA}.fa"
        shell: """
               zcat {input} | fasta2tab | grep "{wildcards.h_lncRNA}" | cut -f2- | tab2fasta -s > {output}
        """
      
rule m_geneName_to_fasta:
        input:
            "/home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M32/transcripts.longest.fa.gz"
        output:
            "best_param_3plex/mm10.{m_lncRNA}.fa"
        shell: """
               zcat {input} | fasta2tab | grep "{wildcards.m_lncRNA}" | cut -f2- | tab2fasta -s > {output}
        """

rule all_fa:
    input:
        lncRNA_human=expand("best_param_3plex/hg38.{h_lncRNA}.fa", h_lncRNA=config["HUMAN_LNCRNAs"]),
        lncRNA_mouse=expand("best_param_3plex/mm10.{m_lncRNA}.fa", m_lncRNA=config["MOUSE_LNCRNAs"])

#change the assemble to select human or mouse, based on the species used
rule bed_to_fasta:
        input:
#            assembly="/home/reference_data/bioinfotree/task/gencode/dataset/hsapiens/32/GRCh38.primary_assembly.genome.clean_id.fa",  #human
            assembly="/home/reference_data/bioinfotree/task/gencode/dataset/mmusculus/M27/GRCm39.primary_assembly.genome.fa", #mouse
            bed_file="{path}.bed"
        output:
                "{path}.bed.fa"
        shell:"""
                bedtools getfasta -name -fi {input.assembly} -bed {input.bed_file} -fo {output}
        """ 

#change the assemble to select human or mouse, based on the species used
    input:
        #lncRNA_human=expand("best_param_3plex/{h_lncRNA}.neg_pos_rand.bed.fa", h_lncRNA=config["HUMAN_LNCRNAs"]),
        lncRNA_mouse=expand("best_param_3plex/{m_lncRNA}.neg_pos_rand.bed.fa", m_lncRNA=config["MOUSE_LNCRNAs"])


#############
# Run 3plex #
#############

REFERENCE_ROOT=os.environ.get("REFERENCE_ROOT")
BIOINFO_REFERENCE_ROOT=REFERENCE_ROOT+"/bioinfotree/task/"
GENCODE_SPECIES=config["GENCODE_SPECIES"]
GENCODE_VERSION=config["GENCODE_VERSION"]
GENCODE_DIR=BIOINFO_REFERENCE_ROOT+"gencode/dataset/"+GENCODE_SPECIES+"/"+GENCODE_VERSION
GENOME_ASSEMBLY=config["GENOME_ASSEMBLY"]
GENOME_ASSEMBLY_DIR=GENCODE_DIR+"/"+GENOME_ASSEMBLY

rule run_3plex:
    input:
        ssRNA="{GeneID}.fa",
        dsDNA="{GeneID}.neg_pos_rand.bed.fa"
    output:
        summary="out_3plex/{GeneID}_ssmasked-{GeneID}.neg_pos_rand.bed.tpx.summary.gz",
	stability="out_3plex/{GeneID}_ssmasked-{GeneID}.neg_pos_rand.bed.tpx.stability.gz"
    shell:"""
        mkdir -p out_3plex; \
        docker run -u `id -u`:`id -g` -it --rm -v $PWD:$PWD imolineris/3plex:v0.1.2-beta $PWD/{input.ssRNA} $PWD/{input.dsDNA} $PWD/out_3plex
    """

rule all_tpx:
    input:
        expand("out_3plex/{GeneID}_ssmasked-dsDNA.tpx.summary.gz", GeneID=config["HUMAN_LNCRNAs"])
    output:
        "ssmasked-dsDNA.tpx.summary_all.gz"
    shell:"""
        matrix_reduce -t 'out_3plex/*_ssmasked-dsDNA.tpx.summary.gz' | gzip > {output}
    """


rule get_coefficients_AUCs:
    output:
        coefficients = "summary_tables/glm_coefficients/{lncRNA}_glm_{stability}-{neg_type}.tsv",
        AUCs = "summary_tables/glm_auc/{lncRNA}_auc_{stability}-{neg_type}.tsv" 
    shell:""" 
        echo 'Please run: get_coefficients+AUCS.Rmd'
    """    

rule coeffiecients_summary:
    input: 
        expand("summary_tables/glm_coefficients/{lncRNA}_glm_{stability}-{neg_type}.tsv", lncRNA=config["LNCRNAs"], neg_type=config["NEG_TYPES"], stability=config["STABILITY"])
    output:
        "summary_tables/coefficients_all.aggregated.tsv"
    shell: """
        matrix_reduce '{input}' -l | grep -v ">" | grep -vw "Intercept" \
        | bawk 'BEGIN {{print "lncRNA", "parameter", "estimate", "p-value", "neg_type", "stability"}} {{print $6, $1, $2, $5, $7, $8}}' > {output} 
    """
    
rule aucs_summary:
    input: 
        expand("summary_tables/glm_auc/{lncRNA}_auc_{stability}-{neg_type}.tsv", lncRNA=config["LNCRNAs"], neg_type=config["NEG_TYPES"], stability=config["STABILITY"])
    output:
        "summary_tables/auc_stab_best.aggregated.tsv"
    shell: """
        matrix_reduce '{input}'| grep -v ">" | grep -vw "Intercept" | grep -v "stability" \
        | bawk 'BEGIN {{print "lncRNA", "AUC_with_shapes", "AUC_stability_only", "AUC_diff", "p-value", "neg_type", "stability"}} {{print$1,$2,$3,$4,$5,$6,$7}}' > {output} 
    """


rule p_adj_tables:
    input:
        ["summary_tables/auc_all.aggregated.tsv", "summary_tables/coefficients_all.aggregated.tsv"]
    output:
        ["summary_tables/auc_all.aggregated.p_adj.tsv", "summary_tables/coefficients_all.aggregated.p_adj.tsv"]
    shell: """
        summary_table/get_summary_tables.R
    """

rule summary_tables_all_AUC:
    input: ['summary_tables/ALL_AUC_stability_best+norm.summary.p-adj.csv', 'summary_tables/ALL_AUC_stability_best.summary.p-adj.csv', 'summary_tables/ALL_AUC_stability_norm.summary.p-adj.csv']
    output: "summary_tables_AUC.all"
    shell: """
        matrix_reduce -t 'summary_tables/*.summary.p-adj.csv' -l '{input}' \
        | bawk 'NR=1{{print}} $2!="lncRNA" {{print}}' > {output}
    """
    
rule summary_tables_all_coeff:
    input: ['summary_tables/ALL_coefficients_no_stability.summary.p-adj.csv', 'summary_tables/ALL_coefficients_stability_best+norm.summary.p-adj.csv', 'summary_tables/ALL_coefficients_stability_best.summary.p-adj.csv', 'summary_tables/ALL_coefficients_stability_norm.summary.p-adj.csv', 'summary_tables/ALL_median_differences.summary.p-adj.csv']
    output: "summary_tables_coeff.all"
    shell: """
        matrix_reduce -t 'summary_tables/*.summary.p-adj.csv' -l '{input}' \
        | bawk 'NR=1{{print}} $2!="lncRNA" {{print}}' > {output}
    """


rule tab2xlsx:
    input: "{file}"
    output: "{file}.xlsx"
    shell: "tab2xlsx < {input} > {output}"


###########################
# Open chromatin fraction #
###########################
rule all_neg_pos_regions:
    input:
        h_reg = expand("best_param_3plex/{h_lnc}.neg_pos_rand.sorted.bed", h_lnc=config["HUMAN_LNCRNAs"]),
        m_reg = expand("best_param_3plex/{m_lnc}.neg_pos_rand.sorted.bed", m_lnc=config["MOUSE_LNCRNAs"])
    output:
        h_pos = "chromatin_fraction/h_pos_regions.bed",
        h_neg = "chromatin_fraction/h_neg_rand_regions.bed",
        m_pos = "chromatin_fraction/m_pos_regions.bed",
        m_neg = "chromatin_fraction/m_neg_rand_regions.bed"
    shell: """
        cat {input.h_reg} | grep "pos" > {output.h_pos};
        cat {input.h_reg} | grep "neg" > {output.h_neg};
        cat {input.m_reg} | grep "pos" > {output.m_pos};
        cat {input.m_reg} | grep "neg" > {output.m_neg}
    """

rule chromatin_fraction:
    input:
        h_pos="chromatin_fraction/h_pos_regions.bed",
        h_neg="chromatin_fraction/h_neg_rand_regions.bed",
        m_pos="chromatin_fraction/m_pos_regions.bed",
        m_neg="chromatin_fraction/m_neg_rand_regions.bed"
    output:
        "chromatin_fraction/chromatin_fraction.txt"
    params:
        h_cCRE="chromatin_fraction/hg38.general_cCRE.bed",
        m_cCRE="chromatin_fraction/mm10.general_cCRE.bed"
    shell: """
        (
        echo "h_pos ({input.h_pos}):" && wc -l {input.h_pos}
        bedtools intersect -loj -a {input.h_pos} -b {params.h_cCRE} | bawk '$7!="."' | cut -f 4 | sort | uniq | wc -l

        echo "h_neg ({input.h_neg}):" && wc -l {input.h_neg}
        bedtools intersect -loj -a {input.h_neg} -b {params.h_cCRE} | bawk '$7!="."' | cut -f 4 | sort | uniq | wc -l

        echo "m_pos ({input.m_pos}):" && wc -l {input.m_pos}
        bedtools intersect -loj -a {input.m_pos} -b {params.m_cCRE} | bawk '$7!="."' | cut -f 4 | sort | uniq | wc -l

        echo "m_neg ({input.m_neg}):" && wc -l {input.m_neg}
        bedtools intersect -loj -a {input.m_neg} -b {params.m_cCRE} | bawk '$7!="."' | cut -f 4 | sort | uniq | wc -l
        ) > {output}
    """
