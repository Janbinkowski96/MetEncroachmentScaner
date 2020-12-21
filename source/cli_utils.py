import click

def print_help(ctx, param, value) -> None:
    """Function to print help msg if no arguments passed into cli"""
    if not value or not param:
        return
    click.echo(ctx.get_help())
    ctx.exit()


def check_flags(ctx, my_norm, sample_sheet, output_directory) -> None:
    """Function to find and print missing parameters"""
    if not my_norm or not sample_sheet or not output_directory:
        arguments = locals().items()
        missing_args = [flag for flag, value in arguments if value is None]
        missing_args = " and ".join(missing_args)
        click.echo(click.style(f"Missing option '{missing_args}', use --help / -h to print options.", bold=True))
        ctx.exit()
        
def break_(msg: str) -> None:
    click.echo(click.style(msg, bold=True))
    sys.exit(-1)


def clear_() -> None:
    click.clear()


def print_(msg: str) -> None:
    click.echo(click.style(msg, bold=True))
