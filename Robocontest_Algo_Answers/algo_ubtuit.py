'''import math

r1 = float(input())
r2 = float(input())
r3 = float(input())

S1 = math.pi * r1**2
S2 = math.pi * r2**2
S3 = math.pi * r3**2

print(f"{S1: .2f}")
print(f"{S2: .2f}")
print(f"{S3: .2f}")

#a = 3.147682823442
#print(f"{a:.4}")

#print(f"{3.1421143541: .2f}")
#import math

#r1 = float(input())
#r2 = float(input())                       #doiraning yuzini topish uchun
#r3 = float(input())

#S1 = math.pi * r1**2
#S2 = math.pi * r2**2
#S3 = math.pi * r3**2

#print(f"{S1:.2f} {S2:.2f} {S3:.2f}")

S = float(input())
h= float(input())                                 #this one "uchburchak yuzini topish "

b = 2*S/h
print(f"{b:.2f}")

import math
r = float(input())
                                                    #this one shar
S =4* math.pi * r *r
print(f"{S: .2f}")
a = int(input())
b = int(input())
c = int(input())                                       #uchburchakning yarimi topish formulasi

P = (a+b+c)/2
print(f"{P:.2f}")

import math
a = float(input())
b = float(input())
h = float(input())

R1 = a/2
R2= b/2                                              # go'laning yuzini topish

l =math.sqrt((R1-R2)**2 + h*h)

S = math.pi*(R1**2 + R2**2 + (R1 +R2)*l)
print(f"{S: .2f}")

import math
h = float(input())
r = float(input())                                    #konus hajmi

V = 1/3 * math.pi * r**2 * h

print(f"{V:.2f}")

v = int(input())
s = int(input())

t = s/v                                                   # masala

print(f"{t:.2f}")
import math

h = float(input())
g = 9.8
t = math.sqrt(2*h/g)                                   #erkin tushush

print(f"{t:.2f}") 1 dan n gacha sonlar berilgan. Berilgan sonlarni yig`indisini toping.

x = int(input())

s = x * 31536000 * (1/1000)                           # 1 sekunda bir milli litr suv

print(f"{s:.0f}")

n = int(input())

S = (n*(n+1))/2

print(f"{S:.0f}")

m = float(input())
g = 9.8
                                                        # fizika
F = m * g
print(f"{F:.2f}")

m = float(input())
a = float(input())

F = m * a                                                 #kuch

print(f"{F:.0f}")
U = int(input())
R = int(input())
                                                           #to'k kuchi
I = U/R

print(f"{I:.3f}")

R1 = float(input())
R2 = float(input())
R3 = float(input())                        #zanjir

R = 1/(1/R1 + 1/R2 + 1/R3)

print(f"{R:.2f}")
import math

x = float(input())
y = float(input())
                                                #chiziqli
c1 = ((x+y)/(y*y + abs((y*y + 2 )/(x + (x**3 / 5) )))) + math.exp(y+2)
print(f"{c1:.2f}")
import math

x = float(input())
y = float(input())                              #chiziqli 5

T11 = ((x**2 + 1)/(x**2+(x*y+y**2)/(y**2+(y+x*y)/(abs(x*y)+5)))) + (1/(1+math.cos(x)+(1/math.sin(abs(x)))))
print(f"{T11:.2f}")
import math
x = float(input())
y = float(input())

f1 = ((2*math.tan(x+(math.pi/6)))/(1/3 + pow(math.cos(y + x **2),1))) + math.log2(x**2+2)
print(f1)
import math

x = float(input())
y = float(input())                              #chiziqli 5

T11 = ((x**2 + 1)/(x**2+(x*y+y**2)/(y**2+(y+x*y)/(abs(x*y)+5)))) + (1/(1+math.cos(x)+(1/math.sin(abs(x)))))
print(f"{T11:.2f}")

import math 
x = float(input())
y = float(input())                             #chiziqli 2 17 misol 

f1 = ((2*math.tan(x+(math.pi/6)))/(1/3 + pow(math.cos(y + x **2),2))) + math.log2(x**2+2)
print(f"{f1:.2f}")

from math import exp
x, y = map(float, input().split())
                                                
c1 = ((x+y)/(y*y + abs((y*y + 2 )/(x + (x**3 / 5))))) + exp(y+2)
print(f"{c1:.2f}")

#from math import exp,cos,atan 
#x = float(input())
#y = float(input()) 

                                                
#f2 = (((1/(x+ (2/x**2) + (3/x**3)))+exp(x**2+3*x))/(atan(x+y)+pow(abs(5+x),2))-pow(cos(y**2+(x**2/2)),2))
#print(f"{f2:.2f}")
#import math 
#x = float(input())

#z = math.log10((x+y)**2+math.sqrt(abs(y)+2)-)
from math import exp,cos,tan,log2,pi
x,y = map(float, input().split())                            

f1 = ((2*tan(x+(pi/6)))/(1/3 + pow(cos(y + x **2),2))) + log2(x**2+2)
print(f"{f1:.2f}")

from math import sin,cos
x,y = map(float, input().split())

T11 = ((x**2 + 1)/(x**2+(x*y+y**2)/(y**2+(y+x*y)/(abs(x*y)+5)))) + (1/(1+cos(x)+(1/sin(abs(x)))))
print(f"{T11:.2f}")

from math import exp,cos,atan 
x,y = map(float, input().split())

                                                
f2 = (((1/(x+ (2/x**2) + (3/x**3)))+exp(x**2+3*x))/(atan(x+y)+pow(abs(5+x),2))-pow(cos(y**2+(x**2/2)),2))
print(f"{f2:.2f}")

from math import log,sqrt,cos,copysign
x,y = map(float, input().split())

z = log(abs(((x+y)**2)+sqrt(abs(y)+2)-(x-(x*y)/((x**2/2)-5))))+ copysign(cos(x+y),2)/pow(x+y,1/3)

print(f"{z:.2f}")
 
a,b = map(float, input().split())

T = (a)**(1/5)+ ((b*(a+b)/(2*b + a*b))**(1/4))*(a**2+b**2+2)        chiziqli 6  21 misol  

print(f"{T:.2f}")  
from math import sin,tan,sqrt,copysign
x1,x2,c,d= map(float, input().split())
c,d = int(c),int(d)
                                                                    chiziqli 7 22 misol
F = abs((sin(abs(c*x2**3+d*x1**3-c*d))**2)/sqrt(c*x1**2+d*x2**2+7))+tan(x1*x2**2+d**3)
print(f"{F:.2f}")
#from math import cos
#a,b,c,d,x= map(float, input().split())
#x = float(x)

#y2 = ((a*x**2+b*x+c)/(x*a**3+a**2+a**(b-c))) + cos(abs((a*x+b)/c*x+d+2**c))

#print(f"{y2:.2f}")

from math import cos,sqrt
a,b,c,x= map(float,input().split())
x = float(x)                             

W2 = 0.75 + ((8.2*x**2+sqrt(abs(x**3+3*x)+cos(x-2))))/((a/4)+(b/3)+(c/2)+1)

print(f"{W2:.2f}")

from math import sqrt,log10
a,x = map(float,input().split())
x = float(x)                                     #chiziqli 10 25 misol
a = int(a)

TT = (sqrt(x-1)+sqrt(x+2)+log10(sqrt(a*x**2)+2))/(sqrt(sqrt(x+2)+sqrt(x+24)+x**5))

print(f"{TT:.2f}")

from math import exp,sin,log,sqrt
a,x,y = map(float,input().split())
a = int(a)
x,y = float(x),float(y)
                                            #chiziqli 11 26 misol 
W2 = sqrt(exp(x*y)-x*sin(a*x)-((x**2+2)/(abs(x)+5)))+sqrt(log(x**2+2)+5)

print(f"{W2:.2f}")

from math import tan,cos,sin,sqrt
x = float(input())
                                        #chiziqli 12 28 misol 
AA = sqrt((2*tan(x+2)-cos(x+2**x))/(1+cos(x+2)**2))+sin(x**2)/(x**2+3)
print(f"{AA:.2f}")
from math import sin,log10,sin,cos
a,x= map(float,input().split())
a = int(a)
x = float(x)                                #chiziqli 13 29 misol

BB1 = x*sin((x/2)+(x/3)+(x/4))+(log10(x**2-2)+3**a)/(cos(x+3)*sin(x+3)+8)

print(f"{BB1:.2f}")

from math import sin,cos,exp,sqrt
x,y,a = map(float,input().split())
a = int(a)  
x = float(x)
y = float(y)                                #chiziqli 14 30 misol hatooooooo

#TT = sqrt(y**2+exp(x)+sqrt(exp(x)+a/(x**2+2))+(cos(x)**2)/(sin(x)**2)+cos(x)**3)
TT = sqrt(y**2 + exp(x) + sqrt(exp(x) + a/(x**2 + 2)) + (cos(x)**2) / (sin(x)**2 + 1e-9) + cos(x)**3)
print(f"{TT:.2f}")
from math import sqrt,sin,exp
x,y,z = map(float,input().split())
x = int(x)
y = float(y)                                  #chiziqli 15 31 misol
z = float(z)

#AF = 2**(-x)*sqrt(x+(abs(y)+2)**1/4)*(((exp(x-1))/(sin(z+2))+2)1/3)
AF = 2**(-x) * sqrt(x + (abs(y) + 2)**(1/4)) * (((exp(x - 1)) / (sin(z + 2)) + 2)**(1/3))
print(f"{AF:.2f}")
from math import cos
a,b,c,d,x= map(float, input().split())
x = float(x)
a = int(a)
b - int(b)                                     chiziqli 8 
c = int(c)
d = int(d)  
y2 = ((a*x**2+b*x+c)/(x*a**3+a**2+a**(b-c))) + cos(abs((a*x+b)/(c*x+d+2**c)))

print(f"{y2:.2f}")

from math import log, sqrt, cos
x , y = map (float, input().split())

z = log(((x+y)**2)+sqrt(abs(y)+2)-(x-((x*y)/((x**2/2)-5)))) + (cos(x+y)**2)/pow(x+y,1/3)
print(f"{z:.2f}")

x, y, z = map(float, input().split())

Max = max(x,y,z)
Min = min(x,y,z)                        #uchta sonning eng katta va eng kichiklarini topish

print(Max,Min)

x,y = map(float,input().split())

Max = max(x,y)
Min = min(x,y)                        #2 sonning eng katta va eng kichiklarini topish

print(Max,Min)
from math import sin,cos,exp,sqrt
x,y,a = map(float,input().split())
a = int(a)  
x = float(x)
y = float(y)                                #chiziqli 14 30 misol hatooooooo

#TT = sqrt(y**2+exp(x)+sqrt(exp(x)+a/(x**2+2))+(cos(x)**2)/(sin(x)**2))+(cos(x)**3)
TT = sqrt(y**2 + exp(x) + sqrt(exp(x) + (a / (x**2 + 2))) + (cos(x)**2 / sin(x**2))) + cos(x)**3
print(f"{TT:.2f}")

from math import cos
a,b,c,d,x= map(float, input().split())
a = int(a)
b - int(b)                                     #chiziqli 8 
c = int(c)
d = int(d) 
if a!=0:
    u  = a**(b-c)
else: u = 1
y2 = ((a*(x**2)+b*x+c)/(x*(a**3)+a**2+u)) + cos(abs((a*x+b)/(c*x+d+2**c)))

print(f"{y2:.2f}")

from math import sin,cos,exp,sqrt
x,y,a = map(float,input().split())
a = int(a)                                 #chiziqli 14 30 misol hatooooooo

#TT = sqrt(y**2+exp(x)+sqrt(exp(x)+a/(x**2+2))+(cos(x)**2)/sin(x**2))+cos(x)**3
TT = sqrt(y**2 + exp(x) + sqrt(exp(x) + a/(x**2 + 2)) + (cos(x)**2) / (sin(x**2))) + cos(x)**3
print(f"{TT:.2f}")

x, y, z = map(float, input().split())

max = x 
min = x
if max<y:
    max = y
if max < z:
    max = z
if min > y:                                ishlad 
    min = y
if min > z:
    min = z
print(max,min)

x, y, z = map(float, input().split()) 
max = x + y + z 
if x>max:
    max = x
if y>max:
    max = y
if z>max:  
    max = z 

min = x+y/2
if x<min:
    min = x
if y<min:
    min = y
if z<min:
    min = z
print(max,min**2)
a ,b ,c = map(int, input().split())
if a<b<c:
    print("YES")                                         #Tarkmoqlanuvchi 4 34 misol 
else:
    print("NO")
a, b, c = map(int, input().split())

if c <= b <= a:
    a*=2
    b*=2
    c*=2

else:                                                    Tarmoqlanuvchi 5 35 misol
    a = abs(a)
    b = abs(b)
    c = abs(c)

print(a, b, c)

a,b = map(int, input().split())
if a>b:                                               Tarqatuvchi 6 36 misol   
    print(a)
else:
    print(a,b)

a ,b = map(int, input().split())

if a <= b:                                           Tarmoqlanuvchi 7 37 misol I should ask from my teacher 
    a = 0

print(a, b)


x, y, z = map(float, input().split())
if 1<=x<=30:
    print(x)

if 1<=y<=30:
    print(y)

if 1<=z<=30:
    print(z)

x,y=map(int,input().split())

if x < y:
    x,y=(x+y)/2,2*x*y
else:                                          wrong answer
    y,x=(x+y)/2,2*x*y

print(x, y)
x,y,z = map(int,input().split())
if x > 0:
    x = x**2
if y > 0:
    y = y**2                                     #Tarmoqlanuvchi 10 40 misol
if z>0:
    z = z**2

print(x,y,z)

x,y,z = map(float,input().split())

if x<1 and y<1 and z<1:                           #tarmoqlanuvchi 11 40 misol 
    min = (x+y)/2
else:
    print(x,y,z)
a,b,c,d = map(float,input().split())
if a<=b<=c<=d:
    max_son = max(a,b,c,d)
    a=b=c=d = max_son                               #tarmoqlanuvchi 12 41 misol a=b=c=d = max_sonfarqii
else:
    min_son = min(a,b,c,d)
    a=b=c=d = min_son

print(a,b,c,d)
x,y=map(float,input().split())

if x<0 and y<0:
    x = abs(x)
    y = abs(y)
    print(x,y)   
if x<0 and y>0 or x>0 and y<0:
    x = x+0.5
    y = y+0.5
    print(x,y)
if x>0 and y>0:
    print(x,y)

x, y = map(float, input().split())

if x < 0 and y < 0:  # Ikkalasi ham manfiy bo'lsa
    x = abs(x)
    y = abs(y)
elif x<0 or y<0:                  #tarmoqlanuvchi 13 43 misol 
    x = x + 0.5
    y = y + 0.5
print(x, y)    

a,b,c = map(int,input().split())

if a+b>c and a+c>b and b+c>a and b+c>a:               #tarmoqlanuvchi 14 44 misol
    print("YES")
else:
    print("NO")
import math

a, b, c = map(int, input().split())
D = b**2 - 4*a*c  

if D >= 0:
    x1 = (-b + math.sqrt(D)) / (2*a)
    x2 = (-b - math.sqrt(D)) / (2*a)             # tarmaoqlanuvchi 15 45 misol
    print(f"{x1:.2f} {x2:.2f}")
else:
    print("NO")
from math import sin
n = int(input())
S = n/2*((sin(1)/2**1)+sin(n)/2**n)
print(f"{S:.2f}")

import math

n = int(input())
S = sum(math.sin(i) / (2**i) for i in range(1, n + 1))
print(f"{S:.2f}")
from math import sin

n = int(input())  
S = 0 

for k in range(1, n + 1):              #sikl 1 62 misol 
    term = sin(k) / (2 ** k) 
    S = S + term 

print(f" {S:.2f}")
from math import sin
n = int(input())                                                   #sikl2 62 msiol 
S = 0 
for i in range(1, n + 1):
    term  = term = ((-1) ** (i-1)) * (sin(i** i)/ (2 ** i))
    S += term
print(f"{S:.2f}")     
from math import factorial
n = int(input())
S = 0 
for i in range(1,n+1):                                                 #sikl3 63 misol
    term = ((-1)**(i-1))*1/(factorial(2*i-1))
    S = S+term
print(f"{S:.4f}") 

n,x= map(float,input().split())
n = int(n)
S = 0
for i in range(1,n+1):                                    #sikl4 64 misol 
    term = ((-1)**(i-1))*(1/(x**(2*i)))
    S += term
print(f"{S:.3f}")

x,n = map(int,input().split())
S = 0 
if x==0:
    print(0)
else:
    for i in range(1,n+1):                       #sikl 65 misol xato
        term = i/(x**(2*i))
        S = S+term
print(f"{S:.3f}")

from math import sin
n,x = map(float,input().split())
n = int(n)                                      #sikl 66 msiol 
S = 0 
for i in range(1,n+1):
    term = ((-1)**(i-1))*(1/i)*sin(i*x)
    S +=term
print(f"{S:.3f}")

from math import sqrt 
n,x = map(float,input().split())
n = int(n)
S = 0                                            #sikl 67 misol
for i in range(1, n+1):
    term = (x**i)/(sqrt(i))
    S+=term
print(f"{S:.3f}")
from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 68 misol
S = 0 
for i in range(1, n+1):
    term = x**i/factorial(i)
    S+=term
print(f"{S:.3f}")
from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 69 misol
S = 0 
for i in range(1, n+1):
    term = ((-1)**i)*((x**i)/factorial(i))
    S+=term
print(f"{S:.3f}")
from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 70 misol
S = 0 
for i in range(1, n+1):
    term = (((-1)**(i-1))*(x**(2*i-1)))/factorial(2*i-1)
    S+=term
print(f"{S:.3f}")

from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 71 misol
S = 0 
for i in range(1, n+1):
    term = (((-1)**(i-1))*(x**(2*i-2)))/factorial(2*i-2)
    S+=term
print(f"{S:.3f}")

from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 72 misol
S = 0 
for i in range(1, n+1):
    term = ((x**(2*i-2)))/factorial(2*i-2)
    S+=term
print(f"{S:.3f}")

n,x = map(float,input().split())
n = int(n)                                        #sikl 73 misol
S = 0 
for i in range(1, n+1):
    term = (x**(2*i-1))/(2*i-1)
    S+=term
print(f"{S:.3f}")

from math import factorial
n,x = map(float,input().split())
n = int(n)                                        #sikl 74 misol
S = 0 
for i in range(1, n+1):
    term = (x**(2*i-1))/(factorial(2*i-1))
    S+=term
print(f"{S:.3f}")

from math import factorial
n,k = map(float,input().split())
n = int(n)                                        #sikl 75 misol  why we give 1 to S 
S = 1 
for i in range(1, n+1):
    term = (((-1)**i)*(k**i))/(factorial(i))
    S+=term
print(f"{S:.3f}")
import math                                          #sikl 76 misol
s = 0
a,b,c = map(int,input().split(' '))
for x in range(a,c+1,3):
    s = s +math.pow((a*x+b)/(b**2+math.cos(x)**2),1/3)-(math.sin(x**2))/(a*b)
print(f"{s:.2f}")

import math
y = 0                                              
a,b,c,d = map(int,input().split())                  #sikl 77 misol
for x in range(c,d+1,2):
    y = y +math.pow((math.sin(a*x)+(b**(2*c)))/(b**2+math.cos(x)**2),1/3) - (math.sin(x**2))/(a*b)
print(f"{y:.2f}")
a,b,c = map(int,input().split())
y = 0                                              #sikl 78 misol
for x in range(a,b+1,2):
    y = y+((a**b)+(b**x)+(c**a))/(2*x**2+3*a**x)
print(f"{y:.2f}")
import math 
a = float(input())
y = 0                                              #sikl 79 misol
for x in range(math.pi/2,math.pi+1,math.pi/19):
    y+= math.pow(a**a,1/3)+x**2*math.cos(a*x)
print(f"{y:.2f}")
    
x , n = map(int,input().split())
S=0
for i in range(n, n+1):
    term = n / (x**2*i)
    S+=term 
print(f"{S:.3f}")
import math 
a = float(input())
y = 0                                              #sikl 79 misol
for x in range(math.pi/2,math.pi+1,math.pi/19):
    y+= math.pow(a**a,1/3)+x**2*math.cos(a*x)
print(f"{y:.2f}")
import math
a = int(input())
h = math.pi/19
x = -math.pi/2
y = 0                                               
while (x<=math.pi):
   y+= math.pow(a**a,1/3)+x**2*math.cos(a*x) 
   x+=h
print(f"{y:.2f}")
a,b,c = map(int,input().split(' '))
x = 5
y = 0
while (x<=10):
   y+=(x**2+b*x+x**c)/(a**2+b**2+x**2)
   x+=0.5                                       # sikl 83 misol 
   print(f"{y:.2f}")
import math 
a = int(input())
h = 0.5 
x = 0                                            # sikl 80 misol 
y = 0 
while (x<=10):
    y += a*math.cos(x)-math.sin(x**2)
    x += h
print(f"{y:.2f}")
import math 
a,b = map(int,input().split( ))
x = 1 
y = 0                                                 #sikl 81 misol 
h = 2 
while (x <= 12):
    y = y + a**2 + ((b + math.sin(x)) / (a**3 + math.cos(x**3)**2)) ** (1/5)
    #y = y + a**2 + (((b+math.sin(x))/(a**3+math.cos(x**3))**2),1/5)
    x += h
print(f"{y:.2f}")

import math
a,b,c = map(int,input( ).split( ))
y = 0 
x = 1                                            #sikl 82 misol 
while x <= 10:
    y +=  (a*x**2)/b + x/c
    x += 3
print(f"{y:.2f}")

a, b, c = map(int,input().split())
h = 0.5
y = 0 
x = 5                                            #sikl 83 misol 

while x <= 10:
    y+= (a**2+b*x+x**c)/(a**2+b**2+x**2)
    x += h 

print(f"{y:.2f}")

import math 
a, b, c = map(int,input().split())
h = 0.25
y = 0 
x = -1                                            #sikl 84 misol 

while x <= 1:
    y+= ((math.sin(a*x)+b**c)/(b**2+math.cos(x)**2))**(1/3) - math.sin(x**2)/(a*b)
    x += h 

print(f"{y:.2f}")


a, b, c = map(int,input().split())
h = 5
y = 0 
x = 1                                            #sikl 85 misol 

while x <= 20:
    y+=(a*x**2+b*x+c)/(a**2+b**2+x**2)
    x += h 

print(f"{y:.2f}")
import math
a, b, c = map(int,input().split())
h = 0.25
x = c 
y = 0
while x <= b:                                      #sikll 86 misol 
    y += a**2*math.cos(x)+ (math.sin(x)/(2)) + b * x**2
    x += h 
print(f"{y:.2f}")
import math 
a = int(input())
h = math.pi/10
y = 0 
x = -math.pi/2                                       # sikl 87 misol
while x <= math.pi :
    y  +=  2 * pow(a**math.sin(2*x),1/3)+ x**2*math.cos(a*x)
    x += h 
print(f"{y:.2f}")

import math 
a, b, c, d = map(int,input().split())
x = d 
h = 1.5 
y = 0                                                #sikl 88 misol
while x <= c:
    y += pow(((a*x+b)/(b**2+math.cos(x)**2)),1/5) - (math.sin(x**2))/(a*b)
    x += h
print(f"{y:.2f}")
import math
a,b,c = map(int,input().split())
x = 0 
h = 0.25 
y = 0                                                 # wring answer sikl 89 misol 
while x <= 1 :                                      
    y += math.sqrt((math.sin(a*x)+b**c)/(b**2+math.cos(x)**2)) - (math.sin(x**2)/(a*b))
    x += h 
print(f"{y:.2f}")

from math import pi, sin, atan, exp, log,log10
a, b, c = map(int,input().split())
x = -pi 
h = pi/10
y = 0                                                #sikl 90 misol
while x <= pi:
    y += (log(pow(a,2*sin(x)))+exp(2*x))/(atan(x)+b)+c
    x += h
print(f"{y:.2f}")

import math 
x,y,a,b = map(int,input().split())
S = 0 
P = 1 
SP = 0 

for k in range(1, x+1):
    S += (k**2+math.sin(k)+5)/((k**7 + 1)**(1/5))
for n in range(1, y+1):
    P*= (n + math.sqrt(n))/(n - (n+1)**(1/5))                           #Sikl 33 
for k in range(1,a+1):
    new = 1 
    for i in range(1,b+1):
        new *= (i**2+(k**2)**(1/i))/((math.sin(i)+math.cos(k))*(i**k))
    SP +=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math
a,b,c,d = map(int,input().split())
S = 0 
P = 1 
SP = 0 

for m in range(1,a+1):
    S += (3*(m**3)+4*m+5)/(m**3+math.log(m))
for k in range(1,b+1):
    P *= k / (k**3+7*k+5)
for i in range(1,c+1):                                            #Sikl 31 misol
    new = 1 
    for m in range(1,d+1):
        new *= (math.log(i)+m**i)/(m**i)
    SP += new
print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math
x,y,a,b = map(int,input().split())
S = 0 
P = 1 
SP = 0 

for a in range(1,x+1):
    S += (a**2 + 2*a)/(a**3 + a* (math.cos(a)**2)+1)
for i in range(1,y+1):
    P*= (i**2 + 1)/(i**(1/3) + 2)
for i in range(1,a+1):                                            #Sikl 32 xatooooo
    new = 1
    for k in range(1,b+1):
        new *= math.log(k**i + k**(1/i)) / (k**3+i**(1/k))
    SP += new
print(f"{S:.2f} {P:.2f} {SP:.2f}")
import math 
x,y,c,d = map(int,input().split())
S = 0
P = 1
SP = 0 

for a in range(1,x+1):
    S += (2*a + math.cos(a))**2 
for a in range(1,y+1):                                       #sikl 34
    P *= (a+6) / math.sqrt(a**2+2)
for k in range(1,c+1):
    new = 0
    for y in range(1,d+1):
        new += (k**2+y) / math.sqrt(k**2+y**2)
    SP +=new 
print(f"{S:.2f} {P:.2f} {SP:.2f}")
import math 
x,y,c,d = map(int,input().split())
S = 0
P = 0
SP = 0 
for i in range(1,x + 1):
    S+= (i**4 + i**2 + 3) / (math.sqrt(i + math.exp(i)))
for k in range(1, y+1):
    P  +=( k +1 )/ (k**3 + 5*k )                                #sikl 35 
for m in range(1, c + 1):
    new = 1 
    for n in range(1, d + 1):
        new *= math.sqrt(abs(m**n - n**m)/ (m **n + n**m) )
    SP +=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")
import math 
x,y,c,d = map(int,input().split())
S = 0
P = 1
SP = 1
for  k in range(1,x + 1):
    S+= ((-1)**k*(k+1))/(k**3 + k **2 +1 )                         #sikl 36
for i in range(1, y + 1):
    P *= (i**3 + abs(i-9)) / (math.log(i)+7*i)
for n in range(1,c+1):
    new = 0
    for m in range(1,d+1):
        new+= ((-1)**m) * ((math.log(m+5))/(m**(n+3)+n*m))
    SP*=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math 
x,y,c,d = map(int,input().split()) 
P = 1 
S = 0 
SP = 0 
for n in range(1, x +1):
    S += 1 / (5 - 17 * n + n**3 )
for m in range(1, y + 1):
    P *= (math.sqrt(abs(m - 5) + 1)) / ( m**2 + 4 * m - 1 )
for i in range(1,c + 1):
    new = 1
    for k in range(1, d + 1):
        new += ((-1)**i) * (((abs(math.sin(k)+math.exp(k)))**(1/7))/(2*abs(4*i**3-k**4)))
    SP+=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math
x,y,c,d = map(int,input().split())
S = 0 
P = 1 
SP = 0 

for a in range(1,x+1):
    S+=(4*a+6*math.log(a))/(a**2+a)
for k in range(1,y+1):
    P*= (abs(a - 6*math.cos(a))) / (a**2 + a**(math.log(a)))
for i in range(1,c+1):                                            
    new = 1 
    for m in range(1,d+1):
        new *=(a*k+x)/(k**2+y**2)
    SP += new
print(f"{S:.2f} {P:.2f} {SP:.2f}")
sonlar = []
sonlarr = []
nnumber = []
number = []
s = 0 
for i in range(5):
    son = int(input(f"{i+1} - soonni kiriting: "))
    sonlar.append(son)  
    s+=son
    a = max(sonlar)
    b = min(sonlar)
    sonlarr.append(son*2)
    if sonlar [i] % 2 == 0:
        nnumber.append(son)
print("eng katta son = ",a,"eng kichik son = ",b)
print(sonlar[::1])
print(len(sonlar))
print(s)
print(sonlarr)
print(nnumber)
import math 
x,y,c,d = map(int,input().split()) 
P = 1 
S = 0 
SP = 0 
for n in range(1, x +1):
    S += 1 / (5 - 17 * n + n**3 )                                # do not use same variable in one 
for m in range(1, y + 1):                                        # 37 misol form sik 
    P *= (math.sqrt(abs(m - 5) + 1)) / ( m**2 + 4 * m - 1 )
for i in range(1,c + 1):
    new = 1
    for k in range(1, d + 1):
        new *= ((-1)**i) * (((abs(math.sin(k)+math.exp(k)))**(1/7))/(2*abs(4*i**3-k**4)))
    SP+=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")
        
import math 
x, y, c, d = map(int, input().split()) 
P = 1 
S = 0 
SP = 0 

for a in range(1, x+1):
    S += (4*a + 6*math.log(a)) / (a**2 + a)

for i in range(1, y+1):                                                  #sikl 38 misol
    P *= abs(i - 6*math.cos(i)) / (i**2 + i**(math.log(i)))

for k in range(1, c+1):
    new = 1
    for b in range(1, d+1):
        new *= (b*k + x) / (k**2 + y**2)
    SP += new

print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math 
x, y, c, d = map(int, input().split()) 
P = 1 
S = 0 
SP = 0 

for k in range(1,x+1):
    S+=k**3+math.exp(k)

for a in range(3,y+1):                                                     #sikl 39 misol 
    P*=(a*x)/(math.sqrt(a**2+x**2))

for i in range(1,c+1):
    new = 1
    for j in range(1,d+1):
        new*=(i*x+j**2)/(math.sqrt(i**2+j*y))
    SP+=new
print(f"{S:.2f} {P:.2f} {SP:.2f}")

import math 
x, y, c, d = map(int, input().split()) 
P = 1 
S = 0 
SP = 1

for a in range(1, x+1):
    S += (a*x + 4) / (math.sqrt(a + math.log(6)))

for b in range(1, y+1):
    P *= (b*x**2 + 6) / (math.sin(b*x))                                #sikl 40 misol 

for i in range(1, c+1):
    for j in range(1, d+1):
        SP*= (i*j + y*x) / (math.sqrt((j*x + y)**i))
    

print(f"{S:.2f} {P:.2f} {SP:.2f}")
n = int(input())
a = list(map(int, input().split()))
s = 0                                                                    #massiv 1 
u = 0
ortacha = sum(a) / n
for i  in range(n):
    if a[i] < ortacha:
        s+=a[i]
        u+=1
print(f"{s/u:.2f}")

n =  int(input())
a_list = list(map(int,input().split()))
a,b = map(int,input().split())
min_elements = min(a_list)

for i in range(a,b+1):
    a_list[i] = a_list[i]/min_elements
    print(a_list[i],end = " ")
'''



        