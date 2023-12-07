### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 6e51732b-4f10-434f-abd1-1af6ca5e3e11
using LinearAlgebra, RowEchelon, Latexify

# ╔═╡ a35ff42e-db2d-47e6-8457-ffa7636677ca
function no_zero_pivots(n)
	A = rand(-9:9, (n,n))
	flag = true
	while flag		
		try 
			lu(A, NoPivot())
		catch 
			A = rand(-9:9, (n,n))
		end
		flag = false
	end
	return convert(Array{Rational}, A)
end;

# ╔═╡ bcc7a290-6817-4746-95a0-7ca22a69c424
function zero_pivots(n)
	flag1 = true
	flag2 = true
	while flag1 || flag2
		A = rand(-9:9, (n,n))
		if det(A) != 0
			flag1 = false
			try 
				lu(A, NoPivot())
			catch
				flag2 = false
				return convert(Array{Rational}, A)
				break
			end
		end
	end
end;

# ╔═╡ 7a54a485-b394-4b18-855f-1b544546f438
md"""
##### Assignment by: Miles Kent
"""

# ╔═╡ a187d952-63a9-11ee-39df-8d663e851841
md"""
# LU Factorization of a Matrix

In this assignment you will need to find the ``LU`` factorization of ``3\times 3`` and ``4\times 4`` matrices. Some with no row exchanges and others which require row exchanges (permutation matrices). 

There is also an extra-credit problem that you should do.

To get credit, you must verify that your factorizations are correct by multiplying them out and showing that their product is the original matrix.

Be sure to **disable the cells** containing the matrices, so their values don't get accidentally changed.
"""

# ╔═╡ 78f18405-2cd0-4cbe-844c-895ffd306c18
md"""
##### 1. Find the ``LU`` factorization of this ``3\times 3`` matrix:
"""

# ╔═╡ 4642d43e-305d-4f78-bf41-844a123dd85c
E31 = [
	1//1 0//1 0//1;
	0//1 1//1 0//1;
	0//1 -2//1 1//1;
]

# ╔═╡ e49fdbd7-571b-4769-8e72-9ce78cc5f734
E21 = [
	1//1 0//1 0//1;
	-2//3 1//1 0//1;
	0//1 0//1 1//1;
]

# ╔═╡ e04e27a8-cafc-49c4-ad90-37e13faad9ad
E32 = [
	1//1 0//1 0//1;
	0//1 1//1 0//1;
	0//1 36//10 1//1;
]

# ╔═╡ b020f1b2-bdf0-488a-9d5b-53bc17e42e09
U = E32*E21*E31*A

# ╔═╡ 981b48a6-3b73-4875-90c9-e58a4b57fac8
L = [
1//1 0//1 0//1;
2//3 1//1 0//1;
4//3 -8//5 1//1
]

# ╔═╡ abbed2b6-1b28-4517-99a8-e66dd83c52e6
L*U==A

# ╔═╡ 132f5ae7-5cec-41fb-9aff-464c0464e726
md"""
##### 2. Find the ``LU`` factorization of this ``3\times 3`` matrix:
"""

# ╔═╡ 18f18956-2a8d-4cf9-abd8-10d9284e5de5
P213 =[
0//1 1//1 0//1;
1//1 0//1 0//1;
0//1 0//1 1//1
]

# ╔═╡ 3f062a18-bed9-4818-b812-95782f51f51e
E31_2 = 
[
1//1 0//1 0//1;
0//1 1//1 0//1;
-7//2 0//1 1//1
]

# ╔═╡ 6d63209c-209e-4d31-b983-3a5715ff4b59
E32_2 = 
[
1//1 0//1 0//1;
0//1 1//1 0//1;
0//1 -27//5 1//1
]

# ╔═╡ 435c1683-948d-4e5f-a542-52dc73fdb763
U_2 = E32_2*E31_2*P213*B

# ╔═╡ a9799c15-dba4-4ce0-ab7b-2a7f8a4c0854
L_2 = [
0//1 1//1 0//1;
1//1 0//1 0//1;
7//2 27//5 1//1
]

# ╔═╡ 150b4656-ecf1-4bdf-b799-4f55910669cd
L_2*U_2==B

# ╔═╡ b0fc4703-dd42-4e29-8ad5-9bc7a9dcdbbd
md"""
##### 3. Find the ``LU`` factorization of this ``4\times 4`` matrix:
"""

# ╔═╡ fc56e1aa-e301-4c84-a9e4-d04c08e3dcc6
E31_3 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
-9//6 0//1 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ dace3490-db71-4756-9ed0-c9902edff26a
E41_3 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
5//6 0//1 0//1 1//1
]

# ╔═╡ 221dc97d-b8f2-469d-be69-a2b4042d2c3f
E32_3 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 -1//2 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ df95aa25-db35-4dcb-8042-80e75b341997
E42_3 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 -1//2 0//1 1//1
]

# ╔═╡ bc899141-2324-4cf4-b06b-0beeda1fb421
E43_3 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 -22//6 1//1
]

# ╔═╡ 1b9da145-67b1-464d-a5a1-6b479cd14b02
U_3 = E43_3*E42_3*E32_3*E41_3*E31_3*C

# ╔═╡ 639f8786-b8c3-44bb-9af3-05bc44e31055
L_3 = [
  1//1  0//1   0//1  0//1;
  0//1  1//1   0//1  0//1;
  3//2  1//2   1//1  0//1;
 -5//6  1//2  11//3  1//1
]

# ╔═╡ 27ec99e5-bb87-40fc-8d5a-f0a61d57b17e
L_3*U_3==C

# ╔═╡ c77a6304-f389-46b0-8fc4-ea96b815cd1c
md"""
##### 4. Find the ``LU`` factorization of this ``4\times 4`` matrix:
"""

# ╔═╡ c69dce16-63f2-404f-8f07-338f1e149b12
P4231_4 = 
[
0//1 0//1 0//1 1//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
1//1 0//1 0//1 0//1
]

# ╔═╡ aa4954db-3824-4c43-a1e3-e99ab6ab1de9
E31_4 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
-2//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ c536c49a-4135-4ffb-ab76-01491f38b231
E21_4 = 
[
1//1 0//1 0//1 0//1;
-5//3 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ 1c8861fb-2f73-4e93-a28f-9d520762d6ad
E32_4 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 -(3*16)//2 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ f1172642-6dc5-4808-88b7-4c77f303708f
E43_4 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 9//333 1//1
]

# ╔═╡ ff53a63d-ec94-47ed-9c60-ad130d7ba14c
U_4 = E43_4*E32_4*E21_4*E31_4*P4231_4*D

# ╔═╡ ce3fbba9-d042-486b-8ec8-f7bbce801c8a
L_4 = [
0//1 0//1 -1//37 1//1;
5//3 1//1 0//1 0//1;
2//1 24//1 1//1 0//1;
1//1 0//1 0//1 0//1
	]

# ╔═╡ 763d4091-c5b1-4cf3-8ca6-4b1b1f2aefe3
L_4*U_4==D

# ╔═╡ 5080dcf5-992e-4d3b-8409-c9d0f32289cc
md"""
##### **Extra Credit**: Find the ``LDU`` factorization of this ``4\times 4`` matrix, where:
* ``L`` is a lower triangular matrix with ``1``'s on the main diagonal.
* ``D`` is a *diagonal matrix* with numbers on the main diagonal and ``0``'s elsewhere.
* ``U`` is an upper triangular matrix with ``1``'s on the main diagonal.
"""

# ╔═╡ 9612e16d-9cc0-45cf-bf77-725a87d2b531
E41_5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 3//1 0//1 1//1
]

# ╔═╡ aee927ef-2c36-4f78-abfc-02c90a9e9d42
E31_5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 3//1 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ cb643ce6-7fa9-4b4d-accf-60fae521fc97
E21_5 = 
[
1//1 0//1 0//1 0//1;
1//8 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ f24f6faf-3961-4ba9-aa23-0f0f1f7b9930
E32_5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 -(8*19)//49 1//1 0//1;
0//1 0//1 0//1 1//1
]

# ╔═╡ 09b3ab47-98b1-4421-b576-05aece3867c0
E42_5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 -(8*15)//49 0//1 1//1
]

# ╔═╡ dcd34a99-7103-4291-88f4-a92df5b34d32
E43_5 = 
[
1//1 0//1 0//1 0//1;
0//1 1//1 0//1 0//1;
0//1 0//1 1//1 0//1;
0//1 0//1 -456//313 1//1
]

# ╔═╡ 81c08f83-24f7-4673-8d95-40af4930a342
U_5 = E43_5*E42_5*E32_5*E21_5*E31_5*E41_5*M

# ╔═╡ f24bf1e1-c131-4d5b-80dd-97d29efb642e
D_5=
	[
8//1 0//1 0//1 0//1;
0//1 -49//8 0//1 0//1;
0//1 0//1 -313//49 0//1;
0//1 0//1 0//1 1779//313
	]

# ╔═╡ 4e23b830-5f9b-4aef-ac80-8744aa6eb465
U_5v2 = [
1//1 -1//8 1//8 7//8
0//1 1//1 -1//49 -15//49
0//1 0//1 1//1 -58//313
0//1 0//1 0//1 1//1
]

# ╔═╡ d97d45be-fbb5-4a77-be7e-cbd4aadfd8f7
U_5 == D_5*U_5v2

# ╔═╡ 490de2bd-64a2-48a4-9285-387beed892bc
L_5 = [
  1//1    0//1     0//1    0//1;
 -1//8    1//1     0//1    0//1;
  3//8    5//49    1//1    0//1;
  3//8  -27//49  456//313  1//1
]

# ╔═╡ 7c416f71-4a9b-446c-b949-4e69f2a7b278
L_5*D_5*U_5v2==M

# ╔═╡ be0404c0-54ce-47f2-b464-586822e970ff
# ╠═╡ disabled = true
#=╠═╡
M = no_zero_pivots(4)
  ╠═╡ =#

# ╔═╡ 83d562e6-8f52-423e-a9f5-c6b7054a6079
M =
[
  8//1  -1//1   1//1  7//1;
 -1//1  -6//1   0//1  1//1;
  3//1  -1//1  -6//1  4//1;
  3//1   3//1  -9//1  9//1	
]

# ╔═╡ 4280220d-8457-4c1a-b458-daf997eedc95
# ╠═╡ disabled = true
#=╠═╡
B = zero_pivots(3)
  ╠═╡ =#

# ╔═╡ 3e997a7f-792f-4aad-9437-bd8559294b3b
C = 
[
 -6//1   0//1  -6//1  -5//1;
  0//1  -4//1   3//1  -9//1;
 -9//1  -2//1  -9//1   1//1;
  5//1  -2//1   1//1   5//1
]

# ╔═╡ fcd2332d-1a3d-4b69-9510-90f57fdf96a3
# ╠═╡ disabled = true
#=╠═╡
A = no_zero_pivots(3)
  ╠═╡ =#

# ╔═╡ 6adab8c9-19c4-4946-96fb-e7642ce398bf
# ╠═╡ disabled = true
#=╠═╡
C = no_zero_pivots(4)
  ╠═╡ =#

# ╔═╡ c0d0bfbc-e231-4202-aa07-033b27ecda3c
B = [
	0//1   5//1  -2//1;
 -2//1  -8//1   9//1;
 -7//1  -1//1   6//1
]

# ╔═╡ 38dc45c8-ea9b-4a4b-99e5-4ceae5e39bca
A = [
-3//1  4//1   5//1;
 -2//1  6//1  -6//1;
 -4//1  0//1  -4//1
]

# ╔═╡ c7f64759-f747-4c99-a6bc-e2ffa29d63fc
# ╠═╡ disabled = true
#=╠═╡
D = zero_pivots(4)
  ╠═╡ =#

# ╔═╡ abc6dc87-6b87-4db7-9d14-be8ac8ea241e
D = 
[
 0//1   0//1  -9//1   2//1;
 5//1  -6//1  -4//1   5//1;
 6//1   8//1   9//1   7//1;
 3//1  -4//1   6//1  -7//1
]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"

[compat]
Latexify = "~0.16.1"
RowEchelon = "~0.2.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "19cc70172396ceef58184fe4a17f4580912aeedf"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "f428ae552340899a935973270b8d98e5a31c49fe"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.1"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RowEchelon]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f479526c4f6efcbf01e7a8f4223d62cfe801c974"
uuid = "af85af4c-bcd5-5d23-b03a-a909639aa875"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"
"""

# ╔═╡ Cell order:
# ╟─6e51732b-4f10-434f-abd1-1af6ca5e3e11
# ╟─a35ff42e-db2d-47e6-8457-ffa7636677ca
# ╟─bcc7a290-6817-4746-95a0-7ca22a69c424
# ╟─7a54a485-b394-4b18-855f-1b544546f438
# ╟─a187d952-63a9-11ee-39df-8d663e851841
# ╟─78f18405-2cd0-4cbe-844c-895ffd306c18
# ╠═fcd2332d-1a3d-4b69-9510-90f57fdf96a3
# ╠═38dc45c8-ea9b-4a4b-99e5-4ceae5e39bca
# ╠═4642d43e-305d-4f78-bf41-844a123dd85c
# ╠═e49fdbd7-571b-4769-8e72-9ce78cc5f734
# ╠═e04e27a8-cafc-49c4-ad90-37e13faad9ad
# ╠═b020f1b2-bdf0-488a-9d5b-53bc17e42e09
# ╠═981b48a6-3b73-4875-90c9-e58a4b57fac8
# ╠═abbed2b6-1b28-4517-99a8-e66dd83c52e6
# ╟─132f5ae7-5cec-41fb-9aff-464c0464e726
# ╠═4280220d-8457-4c1a-b458-daf997eedc95
# ╠═c0d0bfbc-e231-4202-aa07-033b27ecda3c
# ╠═18f18956-2a8d-4cf9-abd8-10d9284e5de5
# ╠═3f062a18-bed9-4818-b812-95782f51f51e
# ╠═6d63209c-209e-4d31-b983-3a5715ff4b59
# ╠═435c1683-948d-4e5f-a542-52dc73fdb763
# ╠═a9799c15-dba4-4ce0-ab7b-2a7f8a4c0854
# ╠═150b4656-ecf1-4bdf-b799-4f55910669cd
# ╟─b0fc4703-dd42-4e29-8ad5-9bc7a9dcdbbd
# ╠═6adab8c9-19c4-4946-96fb-e7642ce398bf
# ╠═3e997a7f-792f-4aad-9437-bd8559294b3b
# ╠═fc56e1aa-e301-4c84-a9e4-d04c08e3dcc6
# ╠═dace3490-db71-4756-9ed0-c9902edff26a
# ╠═221dc97d-b8f2-469d-be69-a2b4042d2c3f
# ╠═df95aa25-db35-4dcb-8042-80e75b341997
# ╠═bc899141-2324-4cf4-b06b-0beeda1fb421
# ╠═1b9da145-67b1-464d-a5a1-6b479cd14b02
# ╠═639f8786-b8c3-44bb-9af3-05bc44e31055
# ╠═27ec99e5-bb87-40fc-8d5a-f0a61d57b17e
# ╟─c77a6304-f389-46b0-8fc4-ea96b815cd1c
# ╠═c7f64759-f747-4c99-a6bc-e2ffa29d63fc
# ╠═abc6dc87-6b87-4db7-9d14-be8ac8ea241e
# ╠═c69dce16-63f2-404f-8f07-338f1e149b12
# ╠═aa4954db-3824-4c43-a1e3-e99ab6ab1de9
# ╠═c536c49a-4135-4ffb-ab76-01491f38b231
# ╠═1c8861fb-2f73-4e93-a28f-9d520762d6ad
# ╠═f1172642-6dc5-4808-88b7-4c77f303708f
# ╠═ff53a63d-ec94-47ed-9c60-ad130d7ba14c
# ╠═ce3fbba9-d042-486b-8ec8-f7bbce801c8a
# ╠═763d4091-c5b1-4cf3-8ca6-4b1b1f2aefe3
# ╟─5080dcf5-992e-4d3b-8409-c9d0f32289cc
# ╠═be0404c0-54ce-47f2-b464-586822e970ff
# ╠═83d562e6-8f52-423e-a9f5-c6b7054a6079
# ╠═9612e16d-9cc0-45cf-bf77-725a87d2b531
# ╠═aee927ef-2c36-4f78-abfc-02c90a9e9d42
# ╠═cb643ce6-7fa9-4b4d-accf-60fae521fc97
# ╠═f24f6faf-3961-4ba9-aa23-0f0f1f7b9930
# ╠═09b3ab47-98b1-4421-b576-05aece3867c0
# ╠═dcd34a99-7103-4291-88f4-a92df5b34d32
# ╠═81c08f83-24f7-4673-8d95-40af4930a342
# ╠═f24bf1e1-c131-4d5b-80dd-97d29efb642e
# ╠═4e23b830-5f9b-4aef-ac80-8744aa6eb465
# ╠═d97d45be-fbb5-4a77-be7e-cbd4aadfd8f7
# ╠═490de2bd-64a2-48a4-9285-387beed892bc
# ╠═7c416f71-4a9b-446c-b949-4e69f2a7b278
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
