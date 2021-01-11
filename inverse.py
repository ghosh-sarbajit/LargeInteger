a=0x89a2359acdb123abdee977
b=0x4589073677236889
inva=a #step 1
invb=b
p=0x7fffffffffffffffffffffffffffffff # 2**127 -1 in hex
# binary value of p-2 '0b1111111111111111111111111111111111111111111111111111
# 111111111111111111111111111111111111111111111111111111111111111111111111101'
for i in range(124):
    inva=(inva*inva)%p
    inva=((inva*a)%p)
    invb=(invb*invb)%p
    invb=((invb*b)%p)

inva=((inva*inva)%p) # for 126 th bit
inva=((inva*inva*a)%p) # for 127 th bit

invb=((invb*invb)%p) # for 126 th bit
invb=((invb*invb*b)%p) # for 127 th bit

print(hex(inva))
print(hex(invb))
