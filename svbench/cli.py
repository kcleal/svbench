import click
import svbench as svb
from svbench import CallSet
import pkg_resources

version = pkg_resources.require("svbench")[0].version


@click.command()
@click.argument('reference_vcf', required=True, type=click.Path())
@click.argument('query_vcfs', required=True, type=click.Path(), nargs=-1)
@click.option("--include", help="Include regions .bed file", type=click.Path(), show_default=False, required=False)
@click.option("--pass-only", help="Assess only PASS variants", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--slop", help="Add intervals +/- slop around breakpoints", default=250, type=int, show_default=True)
@click.option("--min-size-ref", help="Min SV length", default=0, type=int, show_default=True)
@click.option("--min-size-query", help="Min SV length", default=30, type=int, show_default=True)
@click.option("--no-duplicates", help="Don't quantify duplicate true positives", is_flag=True, flag_value=True, show_default=False, default=False)
@click.version_option()
def main(reference_vcf, query_vcfs, include, pass_only, slop, min_size_ref, min_size_query, no_duplicates):
    keep = [svb.Col("FILTER", op="eq", thresh=None)] if pass_only else []
    ref = CallSet(dataset="REFERENCE", no_translocations=False).\
        load_vcf(reference_vcf, other_cols=["FILTER"], keep=keep). \
        filter_by_size(min_size_ref, None, keep_translocations=True)

    if include:
        ref.filter_include_bed(include, inplace=True)
    # print(ref.breaks_df.head())
    ref.add_intervals(slop)
    query = [CallSet(dataset="REFERENCE", caller=path.split("/")[-1], no_translocations=False).
             load_vcf(path, other_cols=["FILTER"], keep=keep).filter_by_size(min_size_query, None, keep_translocations=True)
             for path in query_vcfs]

    if include:
        [i.filter_include_bed(include, inplace=True) for i in query]
    svb.score(ref, query, allow_duplicate_tp=not no_duplicates)


if __name__ == '__main__':
    main()
