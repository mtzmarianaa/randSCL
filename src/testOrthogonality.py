from src.hrhtTest import *

# SCRIPT FOR TESTING THE LOSS OF ORTHOGONALITY IN \hat{H}Vß


testType = 1
#n = 2**7

# Size 1
m = 2**11
n = m
V = createTestMatrices(m, n, type = testType)

err0fro_1, err0_1, sig0_1 = testOrthogonality(V, p = 1)
err1fro_1, err1_1, sig1_1 = testOrthogonality(V, p = 2)
err2fro_1, err2_1, sig2_1 = testOrthogonality(V, p = 2**2)
err3fro_1, err3_1, sig3_1 = testOrthogonality(V, p = 2**3)
err4fro_1, err4_1, sig4_1 = testOrthogonality(V, p = 2**4)
err5fro_1, err5_1, sig5_1 = testOrthogonality(V, p = 2**5)
err6fro_1, err6_1, sig6_1 = testOrthogonality(V, p = 2**6)
err7fro_1, err7_1, _ = testOrthogonality(V, p = 2**7)
err8fro_1, err8_1, _ = testOrthogonality(V, p = 2**8)
err9fro_1, err9_1, _ = testOrthogonality(V, p = 2**9)
err10fro_1, err10_1, _ = testOrthogonality(V, p = 2**10)
err11fro_1, err11_1, _ = testOrthogonality(V, p = 2**11)

errsFRO_1 = [err0fro_1, err1fro_1, err2fro_1, err3fro_1, err4fro_1,
             err5fro_1, err6fro_1, err7fro_1, err8fro_1, err9fro_1, 
             err10fro_1, err11fro_1]
errsSP_1 = [err0_1, err1_1, err2_1, err3_1, err4_1, 
            err5_1, err6_1, err7_1, err8_1, err9_1,
            err10_1, err11_1]

print("Finished size 1")

# Size 2
m = 2**12
n = m
V = createTestMatrices(m, n, type = testType)

err0fro_2, err0_2, sig0_2 = testOrthogonality(V, p = 1)
err1fro_2, err1_2, sig1_2 = testOrthogonality(V, p = 2)
err2fro_2, err2_2, sig2_2 = testOrthogonality(V, p = 2**2)
err3fro_2, err3_2, sig3_2 = testOrthogonality(V, p = 2**3)
err4fro_2, err4_2, sig4_2 = testOrthogonality(V, p = 2**4)
err5fro_2, err5_2, sig5_2 = testOrthogonality(V, p = 2**5)
err6fro_2, err6_2, sig6_2 = testOrthogonality(V, p = 2**6)
err7fro_2, err7_2, _ = testOrthogonality(V, p = 2**7)
err8fro_2, err8_2, _ = testOrthogonality(V, p = 2**8)
err9fro_2, err9_2, _ = testOrthogonality(V, p = 2**9)
err10fro_2, err10_2, _ = testOrthogonality(V, p = 2**10)
err11fro_2, err11_2, _ = testOrthogonality(V, p = 2**11)
err12fro_2, err12_2, _ = testOrthogonality(V, p = 2**12)

errsFRO_2 = [err0fro_2, err1fro_2, err2fro_2, err3fro_2, err4fro_2,
             err5fro_2, err6fro_2, err7fro_2, err8fro_2, err9fro_2, 
             err10fro_2, err11fro_2, err12fro_2]
errsSP_2 = [err0_2, err1_2, err2_2, err3_2, err4_2, 
            err5_2, err6_2, err7_2, err8_2, err9_2,
            err10_2, err11_2, err12_2]

print("Finished size 2")

# # Size 3
# m = 2**13
# n = m
# V = createTestMatrices(m, n, type = testType)

# err0fro_3, err0_3, sig0_3 = testOrthogonality(V, p = 1)
# err1fro_3, err1_3, sig1_3 = testOrthogonality(V, p = 2)
# err2fro_3, err2_3, sig2_3 = testOrthogonality(V, p = 2**2)
# err3fro_3, err3_3, sig3_3 = testOrthogonality(V, p = 2**3)
# err4fro_3, err4_3, sig4_3 = testOrthogonality(V, p = 2**4)
# err5fro_3, err5_3, sig5_3 = testOrthogonality(V, p = 2**5)
# err6fro_3, err6_3, sig6_3 = testOrthogonality(V, p = 2**6)
# err7fro_3, err7_3, _ = testOrthogonality(V, p = 2**7)
# err8fro_3, err8_3, _ = testOrthogonality(V, p = 2**8)
# err9fro_3, err9_3, _ = testOrthogonality(V, p = 2**9)
# err10fro_3, err10_3, _ = testOrthogonality(V, p = 2**10)
# err11fro_3, err11_3, _ = testOrthogonality(V, p = 2**11)
# err12fro_3, err12_3, _ = testOrthogonality(V, p = 2**12)
# err13fro_3, err13_3, _ = testOrthogonality(V, p = 2**13)

# errsFRO_3 = [err0fro_3, err1fro_3, err2fro_3, err3fro_3, err4fro_3,
#              err5fro_3, err6fro_3, err7fro_3, err8fro_3, err9fro_3, 
#              err10fro_3, err11fro_3, err12fro_3, err13fro_3]
# errsSP_3 = [err0_3, err1_3, err2_3, err3_3, err4_3, 
#             err5_3, err6_3, err7_3, err8_3, err9_3,
#             err10_3, err11_3, err12_3, err13_3]

# # Size 4
# m = 2**14
# n = m
# V = createTestMatrices(m, n, type = testType)

# err0fro_4, err0_4, sig0_4 = testOrthogonality(V, p = 1)
# err1fro_4, err1_4, sig1_4 = testOrthogonality(V, p = 2)
# err2fro_4, err2_4, sig2_4 = testOrthogonality(V, p = 2**2)
# err3fro_4, err3_4, sig3_4 = testOrthogonality(V, p = 2**3)
# err4fro_4, err4_4, sig4_4 = testOrthogonality(V, p = 2**4)
# err5fro_4, err5_4, sig5_4 = testOrthogonality(V, p = 2**5)
# err6fro_4, err6_4, sig6_4 = testOrthogonality(V, p = 2**6)
# err7fro_4, err7_4, _ = testOrthogonality(V, p = 2**7)
# err8fro_4, err8_4, _ = testOrthogonality(V, p = 2**8)
# err9fro_4, err9_4, _ = testOrthogonality(V, p = 2**9)
# err10fro_4, err10_4, _ = testOrthogonality(V, p = 2**10)
# err11fro_4, err11_4, _ = testOrthogonality(V, p = 2**11)
# err12fro_4, err12_4, _ = testOrthogonality(V, p = 2**12)
# err13fro_4, err13_4, _ = testOrthogonality(V, p = 2**13)
# err14fro_4, err14_4, _ = testOrthogonality(V, p = 2**14)

# errsFRO_4 = [err0fro_4, err1fro_4, err2fro_4, err3fro_4, err4fro_4,
#              err5fro_4, err6fro_4, err7fro_4, err8fro_4, err9fro_4, 
#              err10fro_4, err11fro_4, err12fro_4, err13fro_4, err14fro_4]
# errsSP_4 = [err0_4, err1_4, err2_4, err3_4, err4_4, 
#             err5_4, err6_4, err7_4, err8_4, err9_4,
#             err10_4, err11_4, err12_4, err13_4, err14_4]

# # Size 5
# m = 2**15
# n = m
# V = createTestMatrices(m, n, type = testType)

# err0fro_5, err0_5, sig0_5 = testOrthogonality(V, p = 1)
# err1fro_5, err1_5, sig1_5 = testOrthogonality(V, p = 2)
# err2fro_5, err2_5, sig2_5 = testOrthogonality(V, p = 2**2)
# err3fro_5, err3_5, sig3_5 = testOrthogonality(V, p = 2**3)
# err4fro_5, err4_5, sig4_5 = testOrthogonality(V, p = 2**4)
# err5fro_5, err5_5, sig5_5 = testOrthogonality(V, p = 2**5)
# err6fro_5, err6_5, sig6_5 = testOrthogonality(V, p = 2**6)
# err7fro_5, err7_5, _ = testOrthogonality(V, p = 2**7)
# err8fro_5, err8_5, _ = testOrthogonality(V, p = 2**8)
# err9fro_5, err9_5, _ = testOrthogonality(V, p = 2**9)
# err10fro_5, err10_5, _ = testOrthogonality(V, p = 2**10)
# err11fro_5, err11_5, _ = testOrthogonality(V, p = 2**11)
# err12fro_5, err12_5, _ = testOrthogonality(V, p = 2**12)
# err13fro_5, err13_5, _ = testOrthogonality(V, p = 2**13)
# err14fro_5, err14_5, _ = testOrthogonality(V, p = 2**14)
# err15fro_5, err15_5, _ = testOrthogonality(V, p = 2**15)

# errsFRO_5 = [err0fro_5, err1fro_5, err2fro_5, err3fro_5, err4fro_5,
#              err5fro_5, err6fro_5, err7fro_5, err8fro_5, err9fro_5, 
#              err10fro_5, err11fro_5, err12fro_5, err13fro_5,
#              err14fro_5, err14fro_5]
# errsSP_5 = [err0_5, err1_5, err2_5, err3_5, err4_5, 
#             err5_5, err6_5, err7_5, err8_5, err9_5,
#             err10_5, err11_5, err12_5, err13_5,
#             err14_5, err15_5]



# Plot how the decrease changes
iter = [1, 2, 2**2, 2**3, 2**4, 2**5, 2**6, 
        2**7, 2**8, 2**9, 2**10, 2**11, 2**12,
        2**13, 2**14, 2**15]
# Relative Frobenius error
plt.figure()
plt.loglog( iter[:12], errsFRO_1, label = 'm = 2**11', c = "#111441" )
plt.loglog( iter[:13], errsFRO_2, label = 'm = 2**12', c = "#A200FF" )
# plt.loglog( iter[:14], errsFRO_3, label = 'm = 2**13', c = "#FF5050" )
# plt.loglog( iter[:15], errsFRO_4, label = 'm = 2**14', c = "#FF9840" )
# plt.loglog( iter, errsFRO_5, label = 'm = 2**15', c = "#004DFF" )
plt.legend()
plt.xlabel("p")
plt.ylabel("Relative Frobenius error")
plt.title("Loss in orthogonality Fro - block randomized Hadamard")

# Relative 2 error
plt.figure()
plt.loglog( iter[:12], errsSP_1, label = 'm = 2**11', c = "#111441" )
plt.loglog( iter[:13], errsSP_2, label = 'm = 2**12', c = "#A200FF" )
# plt.loglog( iter[:14], errsSP_3, label = 'm = 2**13', c = "#FF5050" )
# plt.loglog( iter[:15], errsSP_4, label = 'm = 2**14', c = "#FF9840" )
# plt.loglog( iter, errsSP_5, label = 'm = 2**15', c = "#004DFF" )
plt.legend()
plt.xlabel("p")
plt.ylabel("Relative l2 error")
plt.title("Loss in orthogonality L2 - block randomized Hadamard")


# Zoom into the important part - just using more than p = 2 processors
# Relative Frobenius error
plt.figure()
plt.loglog( iter[1:12], errsFRO_1[1:], label = 'm = 2**11', c = "#111441" )
plt.loglog( iter[1:13], errsFRO_2[1:], label = 'm = 2**12', c = "#A200FF" )
# plt.loglog( iter[1:14], errsFRO_3[1:], label = 'm = 2**13', c = "#FF5050" )
# plt.loglog( iter[1:15], errsFRO_4[1:], label = 'm = 2**14', c = "#FF9840" )
# plt.loglog( iter[1:], errsFRO_5[1:], label = 'm = 2**15', c = "#004DFF" )
plt.legend()
plt.xlabel("p")
plt.ylabel("Relative Frobenius error")
plt.title("Loss in orthogonality Fro - block randomized Hadamard")

# Relative 2 error
plt.figure()
plt.loglog( iter[1:12], errsSP_1[1:], label = 'm = 2**11', c = "#111441" )
plt.loglog( iter[1:13], errsSP_2[1:], label = 'm = 2**12', c = "#A200FF" )
# plt.loglog( iter[1:14], errsSP_3[1:], label = 'm = 2**13', c = "#FF5050" )
# plt.loglog( iter[1:15], errsSP_4[1:], label = 'm = 2**14', c = "#FF9840" )
# plt.loglog( iter[1:], errsSP_5[1:], label = 'm = 2**15', c = "#004DFF" )
plt.legend()
plt.xlabel("p")
plt.ylabel("Relative l2 error")
plt.title("Loss in orthogonality L2 - block randomized Hadamard")



## Plot the change in the singular values
# p = 1
plt.figure()
plt.plot( range(len(sig0_1)),  sig0_1, label = "m = 2**11", linewidth = 0.3, c = "#111441")
plt.scatter( range(len(sig0_1)),  sig0_1, s =  0.4, color = "#111441")
plt.plot( range(len(sig0_2)),  sig0_2, label = "m = 2**12", linewidth = 0.3, c = "#A200FF")
# plt.scatter( range(len(sig0_2)),  sig0_2, s =  0.4, color = "#A200FF")
# plt.plot( range(len(sig0_3)),  sig0_3, label = "m = 2**13", linewidth = 0.3, c = "#FF5050")
# plt.scatter( range(len(sig0_3)),  sig0_3, s =  0.4, color = "#FF5050")
# plt.plot( range(len(sig0_4)),  sig0_4, label = "m = 2**14", linewidth = 0.3, c = "#FF9840")
# plt.scatter( range(len(sig0_4)),  sig0_4, s =  0.4, color = "#FF9840")
plt.legend()
plt.xlabel("Number of singular value")
plt.ylabel("$\sigma_{i}(\hat{H}V)$")
plt.title("Singular values of \hat{H}V with p=$2^0$")

# p = 2^2
plt.figure()
plt.plot( range(len(sig2_1)),  sig2_1, label = "m = 2**11", linewidth = 0.3, c = "#111441")
plt.scatter( range(len(sig2_1)),  sig2_1, s =  0.4, color = "#111441")
plt.plot( range(len(sig2_2)),  sig2_2, label = "m = 2**12", linewidth = 0.3, c = "#A200FF")
plt.scatter( range(len(sig2_2)),  sig2_2, s =  0.4, color = "#A200FF")
# plt.plot( range(len(sig2_3)),  sig2_3, label = "m = 2**13", linewidth = 0.3, c = "#FF5050")
# plt.scatter( range(len(sig2_3)),  sig2_3, s =  0.4, color = "#FF5050")
# plt.plot( range(len(sig2_4)),  sig2_4, label = "m = 2**14", linewidth = 0.3, c = "#FF9840")
# plt.scatter( range(len(sig2_4)),  sig2_4, s =  0.4, color = "#FF9840")
plt.legend()
plt.xlabel("Number of singular value")
plt.ylabel("$\sigma_{i}(\hat{H}V)$")
plt.title("Singular values of \hat{H}V with p=$2^2$")


# p = 5
plt.figure()
plt.plot( range(len(sig5_1)),  sig5_1, label = "m = 2**11", linewidth = 0.3, c = "#111441")
plt.scatter( range(len(sig5_1)),  sig5_1, s =  0.4, color = "#111441")
plt.plot( range(len(sig5_2)),  sig5_2, label = "m = 2**12", linewidth = 0.3, c = "#A200FF")
plt.scatter( range(len(sig5_2)),  sig5_2, s =  0.4, color = "#A200FF")
# plt.plot( range(len(sig5_3)),  sig5_3, label = "m = 2**13", linewidth = 0.3, c = "#FF5050")
# plt.scatter( range(len(sig5_3)),  sig5_3, s =  0.4, color = "#FF5050")
# plt.plot( range(len(sig5_4)),  sig5_4, label = "m = 2**14", linewidth = 0.3, c = "#FF9840")
# plt.scatter( range(len(sig5_4)),  sig5_4, s =  0.4, color = "#FF9840")
plt.legend()
plt.xlabel("Number of singular value")
plt.ylabel("$\sigma_{i}(\hat{H}V)$")
plt.title("Singular values of \hat{H}V with p=$2^5$")


