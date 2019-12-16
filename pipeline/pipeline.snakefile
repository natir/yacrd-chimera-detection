dataset = ["SRR8494940_ont", "SRR8494911_pb", "h_sapiens_chr1", "d_melanogaster_reads": "d_melanogaster_ref", "c_elegans_pb"]
coverage = [i for i in range (1000, 6000, 1000)]
internals = ["1", "0"]
containments = ["1", "0"]
  

rule all:
    [f"yacrd/{reads_name}.yacrd" for reads_name in dataset],
    [f"yacrd/{reads_name}_cov{cov}_i{i}_c{c}_l{length}_.yacrd"
     for reads_name in dataset
     for cov in coverage
     for i in iternals
     for c in containments
    ],
    
rule minimap2_index:
    input:
        ref = "reference/{name}.fasta"
    output:
        index = "reference/{name}.mmi"
    threads:
        16
    params:
        preset = 
    shell:
        "minimap2 -t {threads} -x {params.preset} -d {output.index} {input.ref}"


rule map2ref:
    input:
        reads = "reads/{reads_name}.fasta"
        index = lambda wlc: f"references/{config['reads2ref'][wlc.reads_name]}.mmi"
    output:
        mapping = "mapping/{reads_name}.paf"
    threads:
        16
    params:
        preset = lambda wlc: config['preset'][wlc.reads_name]
    shell:
        "minimap2 -t {threads} -x map-{params.preset} {input.index} {input.reads} > {output.mapping}"


rule reads2reads:
    input:
        reads = "reads/{reads_name}.fasta"
    output:
        ovl = "overlaps/{reads_name}.paf"
    threads:
        16
    params:
        preset = lambda wlc: config['preset'][wlc.reads_name]
    shell:
        "minimap2 -t {threads} -x ava-{params.preset} {input.reads} {input.reads} > {output.ovl}"


rule fpa:
    input:
        all_ovl = "overlaps/{reads_name}.paf"
    output:
        filter_ovl = "overlaps/{reads_name}_i{i}_c{c}_l{length}.paf"
    run:
        cmd = f"fpa -i {wildcards.all_ovl} -o {wilcards.filter_ovl} drop -l {wilcards.length}"

        if wildcards.i == "1":
            cmd += " -i"
        if wildcards.c == "1":
            cmd += " -c"

        shell(cmd)

rule yacrd_basic:
    input:
        ovl = "overlaps/{reads_name}.paf"
    output:
        chimeric_report = "yacrd/{reads_name}.yacrd"
    shell:
        "yacrd chimeric -i {input.ovl} -o {output.chimeric_report}"


rule yacrd_param:
    input:
        ovl = "overlaps/{reads_name}_i{i}_c{c}_l{length}.paf"
    output:
        chimeric_report = "yacrd/{reads_name}_cov{cov}_i{i}_c{c}_l{length}.yacrd"
    shell:
        "yacrd chimeric -c {wildcards.cov} -i {input.ovl} -o {output.chimeric_report}"
