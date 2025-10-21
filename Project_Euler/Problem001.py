def main(n):
    sum = 0
    for i in range(n):
        if (i%3 == 0):
            sum += i
        elif(i%5 == 0):
            sum += i

    print(sum)

main(1000)
