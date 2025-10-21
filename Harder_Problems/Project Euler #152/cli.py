import click
from Problem152 import main

@click.command()
@click.option("--n", "-n", default=80, type=int, help="Maximum denominator.")
@click.option("--print-solutions/--no-print-solutions", default=False, help="Print all solutions.")
def cli(n, print_solutions):
    """Solve Project Euler #152 for a given n."""
    result = main(n, print_solutions)
    click.echo(f"Result for n={n}: {result}")

if __name__ == "__main__":
    cli()