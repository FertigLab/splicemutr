import os

configfile: "config.yaml"

if not os.path.exists(config["GENOTYPES_DIR"]):
    os.mkdir(config["GENOTYPES_DIR"])

rule all:
    input:
        OUT_FILE_GENOTYPES_JSON=config["GENOTYPES_DIR"]+"/genotype_files.txt",
        GENOTYPES_FILE=config["GENOTYPES_DIR"]+"/genotypes.txt",
        GENOTYPES_FILE_FORMATTED=config["GENOTYPES_DIR"]+"/genotypes_reformatted.txt",
        UNIQUE_MHC_FILE=config["GENOTYPES_DIR"]+"/class_1_HLAS.txt"

rule run_arcasHLA:
    input:
        GENOTYPES_DIR=config["GENOTYPES_DIR"],
        FILENAMES_FILE=config["BAMFILES"],
    output:
        OUT_FILE_GENOTYPES_JSON=config["GENOTYPES_DIR"]+"/genotype_files.txt"
    shell:
        """
        START=1
        END=2
        for (( VAR=$START; VAR<=$END; VAR++ ))
        do
        FILE=$(sed -n ${{VAR}}p {input.FILENAMES_FILE})
        FILE_BASE=$(basename $FILE)
        FILE_DIR={input.GENOTYPES_DIR}/${{FILE_BASE}}_dir
        mkdir -p $FILE_DIR

        # sort bam file
        samtools sort -o ${{FILE}}.sorted $FILE

        arcasHLA extract ${{FILE}}.sorted -o $FILE_DIR -v

        cd $FILE_DIR

        FASTQ1=$(ls *.extracted.1*)
        FASTQ2=$(ls *.extracted.2*)
        arcasHLA genotype $FASTQ1 $FASTQ2 -g A,B,C
        -o $FILE_DIR -v
        done
        cd {input.GENOTYPES_DIR}
        find $PWD -type f -name *genotype.json > genotype_files.txt
        """

rule create_genotypes_file:
    input:
        GENOTYPES_FILES=config["GENOTYPES_DIR"]+"/genotype_files.txt"
    output:
        GENOTYPES_FILE=config["GENOTYPES_DIR"]+"/genotypes.txt"
    shell:
        """
        START=1
        END=2
        for (( VAR=$START; VAR<=$END; VAR++ ))
        do
            JSON_FILE=$(sed -n ${{VAR}}p {input.GENOTYPES_FILES})
            echo $(echo $(basename $JSON_FILE) | sed "s/Aligned.genotype.json//g")"\t"$(jq '.[] | .[]' $JSON_FILE | paste -s -d "," | sed 's/"//g') >> {output.GENOTYPES_FILE}
        done
        """

rule format_genotypes_file:
    input:
        GENOTYPES_FILE=config["GENOTYPES_DIR"]+"/genotypes.txt"
    output:
        GENOTYPES_FILE_FORMATTED=config["GENOTYPES_DIR"]+"/genotypes_reformatted.txt",
        UNIQUE_MHC_FILE=config["GENOTYPES_DIR"]+"/class_1_HLAS.txt"
    run:
        class_1_HLAs = set()
        with open(output.GENOTYPES_FILE_FORMATTED,"w") as outfile:
            with open(input.GENOTYPES_FILE) as infile:
                all_lines = infile.readlines()
                for line in all_lines:
                    line_split = line.split("\t")
                    HLAs = line_split[1].split(",")
                    HLAS_shortened = [":".join(HLA.split(":")[0:2]) for HLA in HLAs]
                    HLAS_reformatted = ["HLA-"+(HLA.replace("*","")).replace(":","-") for HLA in HLAS_shortened]
                    class_1_HLAs.update(set([HLA for HLA in HLAS_reformatted if "HLA-A" in HLA or "HLA-B" in HLA or "HLA-C" in HLA]))
                    HLAS_reformatted = ",".join(HLAS_reformatted)
                    line_reformatted=line_split[0]+"\t"+HLAS_reformatted+"\n"
                    outfile.write(line_reformatted)
        with open(output.UNIQUE_MHC_FILE,"w") as outfile:
            class_1_HLAs_list = list(class_1_HLAs)
            class_1_HLAs_string = "\n".join(class_1_HLAs_list)
            outfile.write(class_1_HLAs_string)