from pathlib import Path
import numpy as np
import subprocess
import polars as pl

from typing import NamedTuple

class HMMBuildSettings(NamedTuple):
    match_to_match: float
    match_to_ins: float          
    ins_extend: float
    del_extend: float
    loop_prob: float
    enter_skip_loop: float
    skip_to_skip: float
    match_emit_correct: float

    def to_cli(self):
        return [
            f"--match-to-match {self.match_to_match}",
            f"--match-to-ins {self.match_to_ins}",
            f"--ins-extend {self.ins_extend}",
            f"--del-extend {self.del_extend}",
            f"--loop-prob {self.loop_prob}",
            f"--enter-skip-loop {self.enter_skip_loop}",
            f"--skip-to-skip {self.skip_to_skip}",
            f"--match-emit-correct {self.match_emit_correct}"
        ]



match_to_match_range = np.round(np.linspace(0.85, 0.95, 11, endpoint=True), 4)
match_to_ins_range = np.round(np.linspace(0.01, 0.1, 10, endpoint=True), 4)
ins_extend_range = np.round(np.linspace(0.01, 0.1, 10, endpoint=True), 4)
del_extend_range = np.round(np.linspace(0.05, 0.15, 11, endpoint=True), 4)
loop_prob_range = np.round(np.linspace(0.5, 0.95, 10, endpoint=True), 4)
enter_skip_loop_range = np.round(np.linspace(0.01, 0.3, 10, endpoint=True), 4)
skip_to_skip_range = np.round(np.linspace(0.5, 0.95, 10, endpoint=True), 4)
match_emit_correct_range = np.round(np.linspace(0.5, 0.95, 10, endpoint=True), 4)


np.random.seed(100)

program = "cargo run --release -q --"
in_file = "test/SVA_ref_core.fa"

hyperparam_table = "test/hyper_param_test/hyperparam_table.txt"
meta_f = open(hyperparam_table, "w")
meta_f.write("run\tmatch_to_match\tmatch_to_ins\tins_extend\tdel_extend\tloop_prob\tenter_skip_loop\tskip_to_skip\tmatch_emit_correct\n")


N = int(1e3)
for i in range(N):
    print(f"\r{i}/{N}", end="")

    settings = HMMBuildSettings(
        match_to_match=np.random.choice(match_to_match_range),
        match_to_ins=np.random.choice(match_to_ins_range),
        ins_extend=np.random.choice(ins_extend_range),
        del_extend=np.random.choice(del_extend_range),
        loop_prob=np.random.choice(loop_prob_range),
        enter_skip_loop=np.random.choice(enter_skip_loop_range),
        skip_to_skip=np.random.choice(skip_to_skip_range),
        match_emit_correct=np.random.choice(match_emit_correct_range)
    )

    outfile = f"test/hyper_param_test/test_{i}.txt"
    
    command = f'{program} {in_file} {" ".join(settings.to_cli())} > {outfile} 2> /dev/null'
    subprocess.run(command, shell=True)
    meta_f.write(f"{i}\t{settings.match_to_match}\t{settings.match_to_ins}\t{settings.ins_extend}\t{settings.del_extend}\t{settings.loop_prob}\t{settings.enter_skip_loop}\t{settings.skip_to_skip}\t{settings.match_emit_correct}\n")
meta_f.close()



dfs = []
for i in range(N):
    outfile = f"test/hyper_param_test/test_{i}.txt"
    df = pl.read_csv(outfile, separator='\t', raise_if_empty=False)
    if df.is_empty():
        continue
    dfs.append(df.with_columns(iteration=pl.lit(i)))

total_df = pl.concat(dfs)
total_df.write_csv('test/hyper_param_test/total_results.txt', separator='\t')
for i in range(N):
    outfile = f"test/hyper_param_test/test_{i}.txt"
    Path(outfile).unlink(missing_ok=True)
