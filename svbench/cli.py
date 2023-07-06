import click
import svbench as svb
from svbench import CallSet
import pkg_resources

version = pkg_resources.require("svbench")[0].version


@click.command()
@click.argument('reference_vcf', required=True, type=click.Path())
@click.argument('query_vcfs', required=True, type=click.Path(), nargs=-1)
@click.option("--pass-only", help="Assess only PASS variants", is_flag=True, flag_value=True, show_default=False, default=False)
@click.option("--slop", help="Add intervals +/- slop around breakpoints", default=250, type=int, show_default=True)
@click.option("--min-size", help="Min SV length", default=20, type=int, show_default=True)
@click.version_option()
def main(reference_vcf, query_vcfs, pass_only, slop, min_size):
    keep = [svb.Col("FILTER", op="eq", thresh=None)] if pass_only else []
    ref = CallSet(dataset="REFERENCE", no_translocations=False).\
        load_vcf(reference_vcf, other_cols=["FILTER"], keep=keep). \
        filter_by_size(min_size, None, keep_translocations=True)
    ref.add_intervals(slop)
    query = [CallSet(dataset="REFERENCE", caller=path.split("/")[-1], no_translocations=False).
             load_vcf(path, other_cols=["FILTER"], keep=keep).
             filter_by_size(min_size, None, keep_translocations=True)
             for path in query_vcfs]
    svb.score(ref, query)


if __name__ == '__main__':
    main()
