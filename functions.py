from math import isqrt, ceil
from collections import defaultdict
from typing import List, Dict, Set, Optional

Point = Optional[tuple[int, int]]


class EllipticCurve:
    def __init__(self, p: int, a: int, b: int):
        if p <= 3:
            raise ValueError("p needs to be a number > 3")

        if p % 2 == 0 or p < 2:
            raise ValueError("p needs to be an odd number")
        for i in range(3, isqrt(p) + 1, 2):
            if p % i == 0:
                raise ValueError("p needs to be a prime number")
        
        self.p = p
        self.a = a % p
        self.b = b % p

        discriminant = (4 * pow(self.a, 3, p) + 27 * pow(self.b, 2, p)) % p
        if discriminant == 0:
            raise ValueError(f"The curve is not smooth: 4*{a}³ + 27*{b}² = 0 mod {p}")
        
        self._points: List[Point] = []
        self._order: Optional[int] = None
        self._prime_factors: Dict[int, int] = {}
    
    def _inv_mod(self, x: int) -> int:
        if x % self.p == 0:
            raise ZeroDivisionError(f"Attempt to invert zero mod {self.p}")
        return pow(x, self.p - 2, self.p)
    
    def is_on_curve(self, P: Point) -> bool:
        if P is None: 
            return True
        x, y = P
        if not (0 <= x < self.p and 0 <= y < self.p):
            return False
        left = (y * y) % self.p
        right = (pow(x, 3, self.p) + self.a * x + self.b) % self.p
        return left == right
    
    @property
    def infinity(self) -> Point:
        return None
    
    def negate(self, P: Point) -> Point:
        if P is None:
            return None
        x, y = P
        return (x, (-y) % self.p)
    
    def add(self, P: Point, Q: Point) -> Point:
        if P is None:
            return Q
        if Q is None:
            return P
        
        x1, y1 = P
        x2, y2 = Q

        if x1 == x2 and (y1 + y2) % self.p == 0:
            return None
        
        if P == Q:
            if y1 == 0: 
                return None
            
            numerator = (3 * x1 * x1 + self.a) % self.p
            denominator = (2 * y1) % self.p
            s = numerator * self._inv_mod(denominator) % self.p
        else:
            numerator = (y2 - y1) % self.p
            denominator = (x2 - x1) % self.p
            s = numerator * self._inv_mod(denominator) % self.p

        x3 = (s * s - x1 - x2) % self.p
        y3 = (s * (x1 - x3) - y1) % self.p
        return (x3, y3)
    
    def scalar_mul(self, k: int, P: Point) -> Point:
        if P is None or k == 0:
            return None

        if k < 0:
            k = -k
            P = self.negate(P)
        
        result: Point = None
        current = P
        
        while k > 0:
            if k & 1:  
                result = self.add(result, current)
            current = self.add(current, current) 
            k >>= 1 
        
        return result
    
    def enumerate_points(self, full: bool = True, limit: int | None = None) -> List[Point]:
        if self._points:
            return self._points.copy()
        
        points: List[Point] = []

        for x in range(self.p):
            rhs = (pow(x, 3, self.p) + self.a * x + self.b) % self.p
            for y in range(self.p):
                if (y * y) % self.p == rhs:
                    points.append((x, y))

        points.append(None)
        if full:
            self._points = points.copy()

        if not full and limit is not None:
            return points[:limit]

        return points
    
    
    
    def order(self) -> int:
        if self._order is not None:
            return self._order
        
        points = self.enumerate_points()
        self._order = len(points)
        return self._order
    
    def point_order_naive(self, P: Point) -> int:
        if P is None:
            return 1
        
        current = P
        order = 1
        
        while True:
            current = self.add(current, P)
            order += 1
            if current is None:
                return order
            if order > self.p * 2:
                break
        return self.order()
    
    def point_order_bsgs(self, P: Point) -> int:
        if P is None:
            return 1
        p = self.p
        n_max = p + 1 + 2 * isqrt(p) 
        m = ceil(isqrt(n_max))
        baby_steps = {}
        current = None
        for j in range(m + 1):
            baby_steps[current] = j
            if j < m:
                current = self.add(current, P)
        mP = current
        Q = self.scalar_mul(p + 1, P)
        two_mP = self.scalar_mul(2, mP)
        found = False
        best_k = 0
        best_j = 0
        for k in range(-m, m + 1):
            k_two_mP = self.scalar_mul(k, two_mP)
            R = self.add(Q, k_two_mP)
            if R in baby_steps:
                j = baby_steps[R]
                best_k, best_j = k, j
                found = True
                break
            neg_R = self.negate(R)
            if neg_R in baby_steps:
                j = baby_steps[neg_R]
                best_k, best_j = k, j
                found = True
                break
        if not found:
            return p + 1
        M = p + 1 + 2 * m * best_k - best_j 
        R_check = self.add(Q, self.scalar_mul(best_k, two_mP))
        if R_check in baby_steps:
            pass
        else:
            M = p + 1 + 2 * m * best_k + best_j
        if M <= 0:
            M = abs(M)
        return self._refine_order_correct(P, M)
    
    def _refine_order_correct(self, P: Point, M: int) -> int:
        if M <= 0:
            return 1
        n = self.order()

        divisors = self._get_divisors(n)
        divisors.sort(reverse=True)
        
        for d in divisors:
            if M % d == 0 and self.scalar_mul(d, P) is None:
                return d
        
        factors = self._factorize(M)
        result = M
        
        for prime, exp in factors.items():
            for _ in range(exp):
                candidate = result // prime
                if candidate > 0 and self.scalar_mul(candidate, P) is None:
                    result = candidate
                else:
                    break
        
        return result
    
    def _get_divisors(self, n: int) -> List[int]:
        divisors = set()
        for d in range(1, isqrt(n) + 1):
            if n % d == 0:
                divisors.add(d)
                divisors.add(n // d)
        return list(divisors)
    
    def point_order(self, P: Point) -> int:
        if P is None:
            return 1
        if self.p < 100:
            return self.point_order_naive(P)
        return self.point_order_bsgs(P)
    
    def _factorize(self, n: int) -> Dict[int, int]:
        if n <= 1:
            return {}
        
        factors = defaultdict(int)
        temp = n

        while temp % 2 == 0:
            factors[2] += 1
            temp //= 2
            
        d = 3
        while d * d <= temp:
            while temp % d == 0:
                factors[d] += 1
                temp //= d
            d += 2
        
        if temp > 1:
            factors[temp] += 1
        
        return dict(factors)
    
    def get_group_order_factors(self) -> Dict[int, int]:
        if not self._prime_factors:
            n = self.order()
            self._prime_factors = self._factorize(n)
        return self._prime_factors.copy()
    
    def generate_subgroup(self, P: Point) -> List[Point]:
        if P is None:
            return [None]
        
        order = self.point_order(P)
        subgroup = []
        current = None
        
        for _ in range(order):
            current = self.add(current, P)
            subgroup.append(current)
        
        return subgroup
    
    def find_prime_order_subgroups(self) -> Dict[int, List[List[Point]]]:
        prime_factors = self.get_group_order_factors()
        points = self.enumerate_points()
        
        result: Dict[int, List[List[Point]]] = {}
        seen_subgroups: Set[frozenset] = set()
        
        for prime in prime_factors:
            result[prime] = []
            
            for P in points:
                if P is None:
                    continue
                
                order_P = self.point_order(P)
                if order_P != prime:
                    continue
                
                subgroup = self.generate_subgroup(P)
                subgroup_set = frozenset(subgroup)
                
                if subgroup_set not in seen_subgroups:
                    result[prime].append(subgroup)
                    seen_subgroups.add(subgroup_set)
        
        return {prime: subs for prime, subs in result.items() if subs}
 
   
