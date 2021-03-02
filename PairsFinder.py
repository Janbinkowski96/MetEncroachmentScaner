import pandas as pd
import click

from source.cli_utils import print_help, check_flags, print_
from source.utils import check_cores_number, create_path
from source.process_paralleliser import ParallelProcess
from source.data_processor import DataProcessor
from source.anomaly_scaner import AnomalyScaner


@click.command()
@click.option("--data-type", "-t", type=click.Choice(["E", "K"], case_sensitive=False),
              help="Type of input data EPIC or 450K.")
@click.option("--input-file", "-i", type=click.Path(exists=True), help="Path to myNorm file.")
@click.option("--out", "-o", type=click.Path(dir_okay=True), help="Path to output directory.")
@click.option("--search-range", "-r", type=int, default=50, help="Set range to search for methylation encroachment.")
@click.option("--processes-number", "-c", type=int, default=1,
              help="Set number of cores [more cores make computations faster but use proportional more memory!].")
@click.option('--help', '-h', is_flag=True, expose_value=False, is_eager=False, callback=print_help,
              help="Print help message")
@click.pass_context
def cli(ctx, data_type, input_file, out, search_range, processes_number) -> None:
    check_flags(ctx, data_type, input_file, out)
    check_cores_number(processes_number)

    data_processor = DataProcessor(manifest_type=data_type)
    data_processor.set_manifest()
    types = data_processor.first_load(input_file)
    data_processor.load_mynorm(mynorm_path=input_file, column_types=types)

    data_processor.select_cpg()
    data_processor.split_per_chromosome()
    scanner = AnomalyScaner(search_range=search_range, mynorm=data_processor.mynorm,
                            mynorm_std=data_processor.mynorm_std,
                            annotations=data_processor.manifest)

    process = ParallelProcess(data_processor.cpg_per_chr, processes_number)
    process.run_processes(scanner.search_for_anomalies)
    selected = pd.concat(process.processed, sort=False)

    print_("Setting annotations from EPIC manifest")
    selected = DataProcessor.annotate(selected, data_processor.manifest)

    file_path = create_path(out)
    selected.to_csv(file_path)


if __name__ == '__main__':
    cli()
