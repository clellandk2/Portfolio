max = 4000000

def main():
    sum = 0
    fib = [1,2]
    while (fib[1] < max):
        sum += fib[1]
        fib1 = fib[0] + 2* fib[1]
        fib2 = fib1 + fib[0] + fib[1]
        fib = [fib1, fib2]
    print(sum)

main()
