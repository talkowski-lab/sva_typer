import polars as pl
import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("hmmered_file", type=argparse.FileType("rb"))
args = parser.parse_args()

df = pl.read_csv(args.hmmered_file, separator="\t").filter(
    pl.col("region").is_in(["hexamer_region", "VNTR_region"])
).with_columns(
    length=pl.col("end")-pl.col("start"),
    region=pl.col("region").str.replace("_region", "")
).pivot(
    values=["start", "end", "length"],
    index="ID",
    columns="region"
)
df.columns = [c.replace("_region", "") for c in df.columns]

df.write_csv(sys.stdout, separator="\t")
