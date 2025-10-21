def main():
    count = 4
    primes = [2,3,5,7]
    num = 9
    while (count < 10001):
        prime = 0
        for i in primes:
            if (num%i == 0):
                num = num + 2
                prime = 1
                continue
        if (prime==0): 
            primes.append(num)
            num = num + 2
            count += 1

    print(primes[-1])

main()
