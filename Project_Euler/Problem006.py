n = 100

def main():
    squareSum = (n/2*(n+1))**2
    sumSquares = 0
    for i in range(n):
        j = i + 1
        sumSquares += j**2
    print(squareSum, sumSquares)
    print(squareSum - sumSquares)

main()
