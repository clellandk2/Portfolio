from Problem152 import main  

def test_main_result():
    print("Testing n = 35...")
    assert main(35, False) == 1
    print("\n\nTesting n = 45...")
    assert main(45, False) == 3
    print("\n\nTesting n = 80...")
    assert main(80, False) == 301
    print("\n\nTesting n = 20...")
    assert main(20, False) == 0

    print("\nAll done.\n") 

if __name__ == "__main__":
    test_main_result()