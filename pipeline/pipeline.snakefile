
dataset = ["SRR8494940_ont", "SRR8494911_pb", "h_sapiens_chr1_ont", "d_melanogaster_reads_ont", "c_elegans_pb"]
coverages = [i for i in range(0, 6)]
internals = ["1", "0"]
containments = ["1", "0"]
lengths = list(range(0, 1000, 100)) + list(range(1000, 11000, 1000))

rule all:
    input:
        [f"yacrd/{reads_name}.yacrd" for reads_name in dataset],
        [f"yacrd/{reads_name}_i{i}_c{c}_l{length}.yacrd"
         for reads_name in dataset
         for i in internals
         for c in containments
         for length in lengths
        ],
        [f"yacrd/{reads_name}_i{i}_c{c}_l{length}_cov{cov}.yacrd"
         for reads_name in dataset
         for cov in coverages
         for i in internals
         for c in containments
         for length in lengths
        ],
        [f"yacrd/{reads_name}_g{config['g_distance'][reads_name]}_i{i}_c{c}_l{length}.yacrd"
         for reads_name in dataset
         for i in internals
         for c in containments
         for length in lengths
        ],
        [f"yacrd/{reads_name}_g{config['g_distance'][reads_name]}_i{i}_c{c}_l{length}_cov{cov}.yacrd"
         for reads_name in dataset
         for cov in coverages
         for i in internals
         for c in containments
         for length in lengths
        ],
        [f"mapping/{reads_name}.{config['preset'][reads_name]}.paf" for reads_name in dataset],
    
rule minimap2_index:
    input:
        ref = "references/{name}.fasta"
    output:
        index = "references/{name}.{preset}.mmi"
    threads:
        16
    shell:
        "minimap2 -t {threads} -x map-{wildcards.preset} -d {output.index} {input.ref}"


rule map2ref:
    input:
        reads = "reads/{reads_name}.fasta",
        index = lambda wlc: f"references/{config['reads2ref'][wlc.reads_name]}.{wlc.preset}.mmi"
    output:
        mapping = "mapping/{reads_name}.{preset}.paf"
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
    wildcard_constraints:
        reads_name = "((?!_g).)+"
    threads:
        16
    params:
        preset = lambda wlc: config['preset'][wlc.reads_name]
    shell:
        "minimap2 -t {threads} -x ava-{params.preset} {input.reads} {input.reads} > {output.ovl}"


rule reads2reads_g:
    input:
        reads = "reads/{reads_name}.fasta"
    output:
        ovl = "overlaps/{reads_name}_g{dist}.paf"
    wildcard_constraints:
        dist = "\d+"
    threads:
        16
    params:
        preset = lambda wlc: config['preset'][wlc.reads_name]
    shell:
        "minimap2 -t {threads} -x ava-{params.preset} -g {wildcards.dist} {input.reads} {input.reads} > {output.ovl}"


rule fpa:
    input:
        all_ovl = "overlaps/{reads_name}.paf"
    output:
        filter_ovl = "overlaps/{reads_name}_i{i}_c{c}_l{length}.paf"
    wildcard_constraints:
         length = "\d+"
    run:
        cmd = f"fpa -i {input.all_ovl} -o {output.filter_ovl} drop -l {wildcards.length}"

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
    wildcard_constraints:
        reads_name = "((?!_cov).)+"
    shell:
        "yacrd chimeric -i {input.ovl} -o {output.chimeric_report}"


rule yacrd_param:
    input:
        ovl = "overlaps/{reads_name}.paf"
    output:
        chimeric_report = "yacrd/{reads_name}_cov{cov}.yacrd"
    shell:
        "yacrd chimeric -c {wildcards.cov} -i {input.ovl} -o {output.chimeric_report}"
