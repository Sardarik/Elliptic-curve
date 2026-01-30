from functions import *

print("Hi! This is the instrument to build an elliptic curve over a Galois Field! Please choose if you want to:\n 1. Build an elliptic curve and calculate its order\n"
      " 2. Calculate a point's order\n 3. Find all the subgroups")
first_choice = int(input())
            
            
if first_choice == 1:
    print("Type in p (p > 3):")
    p = int(input())
    print("Type in the first coefficient a:")
    a = int(input())
    print("Type in the second coefficient b:")
    b = int(input())
    curve = EllipticCurve(p, a, b)
    print("Please type in 0 if you want to see every point of the curve or another number if you want to limit the number of points:")
    ec_limit = int(input())
    print(f"The group order is {curve.order()}")
    curve._points = []
    print("Points of the group:")
    if ec_limit==0: print(curve.enumerate_points())
    else: print(curve.enumerate_points(full=False, limit=3))
    

if first_choice == 2:
    print("Type in p (p > 3):")
    p = int(input())
    print("Type in the first coefficient a:")
    a = int(input())
    print("Type in the second coefficient b:")
    b = int(input())
    curve = EllipticCurve(p, a, b)
    
    print("Type in the point. First enter the x-axis coordinate, then the y-axis coordinate:")
    p_x = int(input())
    p_y = int(input())
    P = (p_x, p_y)

    if not curve.is_on_curve(P):
        print(f"Error: point {P} is not on the curve!")
    else:
        order = curve.point_order(P)
        print(f"Point order: {order}")

        print("\nNow enter the multiplicity k (integer) to compute kP:")
        k = int(input())
        result = curve.scalar_mul(k, P)
        
        if result is None:
            print(f"{k}P = O (point at infinity)")
        else:
            print(f"{k}P = ({result[0]}, {result[1]})")


if first_choice == 3:
    print("Type in p (p > 3):")
    p = int(input())
    print("Type in the first coefficient a:")
    a = int(input())
    print("Type in the second coefficient b:")
    b = int(input())
    curve = EllipticCurve(p, a, b)
    subgroups = curve.find_prime_order_subgroups()
    
    if subgroups:
        for prime, subgr_list in subgroups.items():
            print(f"\nSubgroups of order {prime}:")
            print(f"  Number of such subgroups: {len(subgr_list)}")

            for i, subgroup in enumerate(subgr_list, 1):
                print(f"\n  Subgroup {i}:")
                for j, point in enumerate(subgroup):
                    if point is None:
                        print(f"    {j+1:2d}. O (point at infinity)")
                    else:
                        print(f"    {j+1:2d}. ({point[0]}, {point[1]})")
    else:
        print("\nNo prime order subgroups found!")
 
    
