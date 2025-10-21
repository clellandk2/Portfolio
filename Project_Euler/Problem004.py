    #Intuition says the largest palindrome made from the product of two 3-digit numbers will be 6 digits, 
    #and likely begin with a 9, and therefore likely end with a 9, and therefore, be odd. 
    #So it will be the product of two odd numbers. I'll use that to narrow my search space. 

def factors(n):
    div = 999
    div2 = 1
    while (div > 99): #this should never be the failing condition, but just to guard against infinte loop glitches
        div2 = n / div
        if (div2 > 999):
            break
        if (n%div == 0):
            print(div, div2)
            return 1
        else: 
            div = div -2
    return 0

def main():
    max = 998001  #(999*999)

    for a in range(9,-1,-1):
        for b in range(9,-1,-1):
            n = 900009 + a * 10010 + b * 1100
            if (factors(n) == 1):
                print(n)
                return

main()
