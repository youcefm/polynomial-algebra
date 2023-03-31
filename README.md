# polynomial-algebra
Seamless computation of basic Polynomial operations

## Example Code

```python
from polynomial import Polynomial
poly1 = Polynomial([1,2,3])
poly2 = Polynomial([0,1,0,4])
print(poly1, poly2)
print('evaluate poly1 at x=5: ', poly1(5))
print('evaluate poly2 at x=5: ', poly2(5))
print(poly1 + poly2)
print(poly1 - poly2)
print(poly1 * poly2)
print(poly2 / poly1)
print(poly2.derivative())
print(poly2.find_roots())
```
