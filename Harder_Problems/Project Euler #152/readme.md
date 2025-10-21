# Project Euler #152 – "Sums of Square Reciprocals"

This directory contains a Python solution to Project Euler problem #152.  
By default, the code solves the **n = 80** formulation, but n can easily be changed either via the **Click CLI** or by modifying the pytest file.  
On my laptop, the runtime for n=80 is approximately **2 seconds**.

> Note: For n > 81, additional primes larger than 37 would need to be added manually to the main script.  
> Also, the runtime grows exponentially with n, so proceed with caution. The code does **not** implement a sieve of Eratosthenes for generating primes automatically.

---

## Usage

Run the program using the CLI with **Click**. There are two parameters:

- `--n` or `-n` : maximum denominator value (before squaring).  
- `--print-solutions` / `--no-print-solutions` : if true, prints all valid solutions in addition to the total count.

**Example:**

```bash
python cli.py --n 45 --print-solutions
python cli.py -n 80

## Testing

A pytest file is included to verify correctness:

Checks the example cases n = 35 and n=45
Checks the known solution for n=80
Checks that no solutions are found for n = 20

Run tests via:
python test_problem152.py
# or if pytest is installed
pytest test_problem152.py
