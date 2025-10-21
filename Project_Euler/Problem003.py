def main():
    n = 600851475143
    i = 3
    largest = 1
    while (i*i < n): 
        if (n%i == 0):
            n = n//i
            largest = i
        else:
           i+= 2 
    print(n)

main()
