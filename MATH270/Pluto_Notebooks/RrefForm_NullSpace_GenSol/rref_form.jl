### A Pluto.jl notebook ###
# v0.19.29

using Markdown
using InteractiveUtils

# ╔═╡ 3cdc1efc-8382-4242-bfe4-5a5f5cff36b2
using LinearAlgebra, RowEchelon, Latexify, StatsBase, Shuffle, LaTeXStrings

# ╔═╡ 442e105b-04c0-4b11-8a9c-938a9bb61654
function shuffle_cols(M)
	return (M'[shuffle(1:end), :])'
end;

# ╔═╡ b4ddfea3-ce79-49e7-bd66-4fb631d44944
function id(n)
	return Matrix{Rational{Int64}}(I,n,n)
end;

# ╔═╡ ac69423d-3c0e-4895-bd5c-ed0193f6c059
function randmat1(m,n,r)
	A = rand(-9//1:9//1, (m,n))
	B = zeros(m,r)
	b = A*rand(-9//1:9//1, n)
	b = reshape(b, (length(b), 1))
	
	for k in 1:r
		combnum = rand(1:n)
		cols = sample(1:n, combnum, replace=false)
		# println(cols)
		newcol = zeros(m,1)
		for col in cols
			coeff = 0
			while coeff == 0
				coeff = rand(-5:5)
			end
			newcol += coeff * A[:,col]
		end
		B[:,k] += newcol
		B = [floor(Int, x)//1 for x in B]
	end
	M = [A B]
	M = shuffle_cols(M)
	return M, b
end;

# ╔═╡ 4826f03d-d9b7-405c-997e-aabdc82a2ba8
function randmat()
	m = rand(3:4)
	n = rand(2:4)
	r = rand(2:3)

	randmat1(m,n,r)
end;

# ╔═╡ 8485ccc0-b25e-4842-aa74-d172a83f8915
function testsoln(A,X,xp,b)
	test=[]
	numcols = size(X)[2]
	b = vec(b)
	for k in 1:5
		coeffs = rand(-20//1:20//1, (numcols,1))
		xn = vec(X*coeffs) # vec() converts nx1 matrix to a vector of length n
		xgen = xp + xn # xgen of a vector
		tst = A*xgen == b
		push!(test, tst)
	end
	return test
end;

# ╔═╡ 459b89ea-6933-11ee-3799-0decb576a655
md"""
##### Assignment by: Miles Kent
"""

# ╔═╡ 3a7397dc-7f25-48d9-85c3-2d8c4d79e296
md"""
# Reduced Row Echelon Form

When we reduce an ``n\times m`` matrix ``A`` to reduced row echelon (rref) form, no only do we learn about the null space of ``A``, but we also learn how the dependent columns of ``A`` are linear combinations of the independent (pivot) columns and how to factor ``A`` as ``CR``, where ``C`` contains the independent columns of ``A`` and ``R = \text{rref}(A)``.

In this assignment you can use the function "rref" on a matrix to reduce it to reduced row echelon form. For example, if

```math
A = \begin{bmatrix*}[r]
	2 & -1 & 0 & 3 & 1\\
	-1 & 4 & -3 & 1 & 1\\
	3 & -1 & 2 & -2 & 4
	\end{bmatrix*}
```
then 

```math
\text{rref}(A) = \begin{bmatrix}
	1 & 0 & 0 & \frac{11}{17} & \frac{19}{17} \\
	0 & 1 & 0 & \frac{-29}{17} & \frac{21}{17} \\
	0 & 0 & 1 & \frac{-48}{17} & \frac{16}{17} \\
	\end{bmatrix}
```
"""


# ╔═╡ 12fb7a36-ad71-4772-ac41-e42e9e8cf2d8
md"""
#### Helper Functions
I've created several helper functions to aid in completing this assignment. They are:

* randmat(): returns a matrix ``A`` and a vector ``\pmb{b}`` in the equation ``A\pmb{x} = \pmb{b}`` for your problem.
* id(n): Returns the ``n\times n`` identity matrix
* testsoln(A,X,xp,b): Here ``A`` is the matrix returned by randmat(), ``X`` is the matrix of null space basis vectors, xp is a particular solution to ``A\pmb{x} = \pmb{b}``, and ``\pmb{b}`` is the vector returned by the randmat() function. This function tests five random solution vectors ``\pmb{x} = \pmb{x}_p + \text{(linear combination of nullspace vectors)}`` in the equation ``A\pmb{x} = \pmb{b}``. If the vector solves the equation, it returns "true", otherwise it returns "false".
"""

# ╔═╡ a4368e61-e504-44c0-ae5b-aac5c8294741
md"""
For each problem below, you must:
* Find matrices ``C`` and ``R``, then show that ``A = CR``
* Find a particular solution to ``Ax = b``
* Find the matrix ``X`` of null space basis vectors
* Use the testsoln function to check your answers. You should get five "true" results.
"""

# ╔═╡ 79100d02-6fdc-49a2-8108-fa8aa849e52c
md"""
##### Example

Consider the matrix ``A0`` below:
"""

# ╔═╡ a0833cfb-510c-4864-ab75-322835423dc0
# ╠═╡ disabled = true
#=╠═╡
A0, b0 = randmat();
  ╠═╡ =#

# ╔═╡ 4857ff2d-0e9b-4edf-83d6-1aa5d1b391f0
A0 = [-8//1   8//1  -16//1   32//1;
       5//1  -7//1   14//1  -30//1;
       3//1  -4//1    8//1  -17//1]; latexify(A0)

# ╔═╡ c071ce04-7943-4003-b891-0544fbe1207d
b0 = [-16//1; 20//1; 11//1]; latexify(b0)

# ╔═╡ 9b0cb42b-0700-4103-b89b-1036cf031ee3
A0b0 = [A0 b0]; latexify(A0b0) # Augmented matrix

# ╔═╡ 07fc1314-19f0-4002-b8b2-246f2326c013
A0b0rr = rref(A0b0); latexify(A0b0rr) # Reduced row echelon form of augmented matrix

# ╔═╡ 90c7dd8a-cc98-4edb-8fe3-2462cbdfe10f
A0rr = A0b0rr[:,1:4]; latexify(A0rr) # 

# ╔═╡ 6a11cac8-aacf-4b26-96fc-7215169e3122
md"""
##### The matric ``C``
"""

# ╔═╡ 47f9cbe0-ed42-4eb9-83a3-db193aa464aa
C = A0[:,[1,2]]; latexify(C) # Columns are basis of column space of A0

# ╔═╡ bd816436-6bef-47e2-a0a1-49959a5da653
md"""
##### The matrix R
"""

# ╔═╡ b612b778-854d-4556-879b-e0cf2ee8cb02
Arr = A0rr[1:2,:]; latexify(Arr) # Reduced row echelon form of A0

# ╔═╡ 472c021d-713d-4d84-b613-7c1b96778497
md"""
##### Checking that ``A = CR``
"""

# ╔═╡ 1b86880d-07c5-4817-ace6-96403e5788fd
A0==C*Arr # A0 factored as C * Arr

# ╔═╡ 10f72d54-8676-4af9-8006-ee90a11bfab7
md"""
##### A particular solution ``x_p`` to ``Ax = b``
"""

# ╔═╡ 9b409192-5c4e-4fb4-9f83-4dbf4edfa293
xp = [-3//1; -5//1; 0; 0] # particular solution to Ax = b

# ╔═╡ 6f2574fb-d012-4dc8-bdbc-e507f1d8873e
F = Arr[:,[3,4]]; latexify(F)

# ╔═╡ 40875b43-4acc-4903-a16e-0e02c2cd246d
md"""
##### The matrix of null space basis vectors; each column is a basis vector of ``N(A)``
"""

# ╔═╡ 341d0a84-da4d-4be9-a096-077a9a2386a3
X = [-F; id(2)]; latexify(X) # Columns are basis of null space

# ╔═╡ a33bf853-40a0-44ee-80d7-82c3d3dab5ab
md"""
##### Using the "testsoln" function to check our answers
"""

# ╔═╡ 1f73ea45-ab0e-4ce6-909e-8dc63453a4b0
testsoln(A0,X,xp,b0) # Check that general solution satisfies Ax = b

# ╔═╡ e2226881-26a2-48d4-8432-95363aeee1dc
md"""
### The Problems
For each problem, copy the ``A`` matrices and the ``b`` vectors into the same variables after disabling the cell containing the randmat() function.
"""

# ╔═╡ 5ebcceff-6597-4cd8-ad9b-fc0c68501833
md"""
##### Problem 1
"""

# ╔═╡ 3d1cbd4b-c5c2-4d11-bcb0-bfb9db3d6ee4
Aug1 = [A1 b1]; latexify(Aug1)

# ╔═╡ 21d0e6e2-239f-4d0a-841c-456a1a8c913b
rrefAug1 = rref(Aug1); latexify(rrefAug1)

# ╔═╡ c07f0afb-c31f-4f2a-b44d-107d99e30dd1
md"""
 #### Matrix C
"""

# ╔═╡ 756f1ae9-9e7d-4cac-a642-c91c8ace0294
C1 = A1[:,1:2]; latexify(C1)

# ╔═╡ 9a5b59f9-8cdd-418b-8759-ce26c4ba6d76
md"""
#### Matrix R
"""

# ╔═╡ 8855270c-1e5b-437b-9cf1-4255317bbc94
R1 = rrefAug1[1:2,1:4]; latexify(R1)

# ╔═╡ 4d1acbd4-f5b1-401b-a3fb-d6791ea01ac1
md"""
#### Verify $A = CR$
"""

# ╔═╡ 4ebae229-f9c4-4999-b728-55147253ae04
A1==C1*R1

# ╔═╡ 32eab404-aee5-4e20-acbb-c6510cb205bd
md"""
#### Particular Solution $x_p$
"""

# ╔═╡ cb3d36f9-5e01-4c37-b68a-74817bcc52d8
xp1 = [-3//1; 5//1; 0//1; 0//1]; latexify(xp1)

# ╔═╡ 76d8cd59-9225-4aee-b180-cbdca5bbf11c
A1*xp1==b1

# ╔═╡ bd1c6a0b-a464-40c7-b139-020f29ef6a57
I2_1 = id(2); latexify(I2_1)

# ╔═╡ 74f40b6d-d36a-45ea-8a74-9842873f7480
F1 = R1[:,[3,4]]; latexify(F1)

# ╔═╡ 7c4f0fb5-86b1-4245-bc80-f45e7dfe356f
md"""
Null Space null(A1) basis vectors matrix
"""

# ╔═╡ b8d9d479-4f0e-4ad9-806c-3659966d7bef
md"""
##### Problem 2
"""

# ╔═╡ 2e7b97dc-79fa-4b11-8af3-cc5c21d536d8
A2, b2 = randmat(); latexify(A2)

# ╔═╡ c3cedb6d-27c1-4788-a7d9-1575b0d7c2e9
latexify(b2)

# ╔═╡ dc6e4adb-d7d3-4d35-afe3-7c7eb1f3ec51
md"""
##### Problem 3
"""

# ╔═╡ 83724a2a-b8ec-4008-a02c-e7ff535f6f87
A3, b3 = randmat(); latexify(A3)

# ╔═╡ 2189aa54-5ec3-4da7-abf1-86633bcc3b3d
latexify(b3)

# ╔═╡ d44f7b94-e2e5-4c92-970b-32ffffc81fa1
# ╠═╡ disabled = true
#=╠═╡
A1, b1 = randmat(); 
  ╠═╡ =#

# ╔═╡ 6bbd0cda-2ec7-4dbe-82a6-3cafb7c2ea99
b1 = [  3//1;
 61//1;
 43//1]

# ╔═╡ 097bcd12-233a-40ab-bba4-351ff7398461
A1 = [-6//1  -3//1   6//1  21//1;
 -32//1  -7//1  -4//1  31//1;
 -26//1  -7//1   2//1  37//1]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
RowEchelon = "af85af4c-bcd5-5d23-b03a-a909639aa875"
Shuffle = "bf21e494-c40e-4daa-abfb-de5ec0aad010"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.16.1"
RowEchelon = "~0.2.1"
Shuffle = "~0.1.1"
StatsBase = "~0.34.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.3"
manifest_format = "2.0"
project_hash = "c921dbf506c00d1ba6b2de783f964a6a0d5c58c4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Compat]]
deps = ["UUIDs"]
git-tree-sha1 = "8a62af3e248a8c4bad6b32cbbe663ae02275e32c"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.10.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+0"

[[deps.DataAPI]]
git-tree-sha1 = "8da84edb865b0b5b0100c0666a9bc9a0b71c553c"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.15.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3dbd312d370723b6bb43ba9d02fc36abade4518d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.15"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

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

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.OrderedCollections]]
git-tree-sha1 = "2e73fe17cac3c62ad1aebe70d44c963c3cfdc3e3"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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

[[deps.Shuffle]]
deps = ["Random"]
git-tree-sha1 = "b812fb30d6d8b295b71dd5a4102d1ae7b60698e3"
uuid = "bf21e494-c40e-4daa-abfb-de5ec0aad010"
version = "0.1.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "c60ec5c62180f27efea3ba2908480f8055e17cee"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─3cdc1efc-8382-4242-bfe4-5a5f5cff36b2
# ╟─442e105b-04c0-4b11-8a9c-938a9bb61654
# ╟─b4ddfea3-ce79-49e7-bd66-4fb631d44944
# ╟─4826f03d-d9b7-405c-997e-aabdc82a2ba8
# ╟─ac69423d-3c0e-4895-bd5c-ed0193f6c059
# ╟─8485ccc0-b25e-4842-aa74-d172a83f8915
# ╟─459b89ea-6933-11ee-3799-0decb576a655
# ╟─3a7397dc-7f25-48d9-85c3-2d8c4d79e296
# ╟─12fb7a36-ad71-4772-ac41-e42e9e8cf2d8
# ╟─a4368e61-e504-44c0-ae5b-aac5c8294741
# ╟─79100d02-6fdc-49a2-8108-fa8aa849e52c
# ╠═a0833cfb-510c-4864-ab75-322835423dc0
# ╠═4857ff2d-0e9b-4edf-83d6-1aa5d1b391f0
# ╠═c071ce04-7943-4003-b891-0544fbe1207d
# ╠═9b0cb42b-0700-4103-b89b-1036cf031ee3
# ╠═07fc1314-19f0-4002-b8b2-246f2326c013
# ╠═90c7dd8a-cc98-4edb-8fe3-2462cbdfe10f
# ╟─6a11cac8-aacf-4b26-96fc-7215169e3122
# ╠═47f9cbe0-ed42-4eb9-83a3-db193aa464aa
# ╟─bd816436-6bef-47e2-a0a1-49959a5da653
# ╠═b612b778-854d-4556-879b-e0cf2ee8cb02
# ╟─472c021d-713d-4d84-b613-7c1b96778497
# ╠═1b86880d-07c5-4817-ace6-96403e5788fd
# ╟─10f72d54-8676-4af9-8006-ee90a11bfab7
# ╠═9b409192-5c4e-4fb4-9f83-4dbf4edfa293
# ╠═6f2574fb-d012-4dc8-bdbc-e507f1d8873e
# ╟─40875b43-4acc-4903-a16e-0e02c2cd246d
# ╠═341d0a84-da4d-4be9-a096-077a9a2386a3
# ╟─a33bf853-40a0-44ee-80d7-82c3d3dab5ab
# ╠═1f73ea45-ab0e-4ce6-909e-8dc63453a4b0
# ╟─e2226881-26a2-48d4-8432-95363aeee1dc
# ╟─5ebcceff-6597-4cd8-ad9b-fc0c68501833
# ╠═d44f7b94-e2e5-4c92-970b-32ffffc81fa1
# ╠═097bcd12-233a-40ab-bba4-351ff7398461
# ╠═6bbd0cda-2ec7-4dbe-82a6-3cafb7c2ea99
# ╠═3d1cbd4b-c5c2-4d11-bcb0-bfb9db3d6ee4
# ╠═21d0e6e2-239f-4d0a-841c-456a1a8c913b
# ╟─c07f0afb-c31f-4f2a-b44d-107d99e30dd1
# ╠═756f1ae9-9e7d-4cac-a642-c91c8ace0294
# ╟─9a5b59f9-8cdd-418b-8759-ce26c4ba6d76
# ╠═8855270c-1e5b-437b-9cf1-4255317bbc94
# ╟─4d1acbd4-f5b1-401b-a3fb-d6791ea01ac1
# ╠═4ebae229-f9c4-4999-b728-55147253ae04
# ╟─32eab404-aee5-4e20-acbb-c6510cb205bd
# ╠═cb3d36f9-5e01-4c37-b68a-74817bcc52d8
# ╠═76d8cd59-9225-4aee-b180-cbdca5bbf11c
# ╠═bd1c6a0b-a464-40c7-b139-020f29ef6a57
# ╠═74f40b6d-d36a-45ea-8a74-9842873f7480
# ╟─7c4f0fb5-86b1-4245-bc80-f45e7dfe356f
# ╟─b8d9d479-4f0e-4ad9-806c-3659966d7bef
# ╠═2e7b97dc-79fa-4b11-8af3-cc5c21d536d8
# ╠═c3cedb6d-27c1-4788-a7d9-1575b0d7c2e9
# ╟─dc6e4adb-d7d3-4d35-afe3-7c7eb1f3ec51
# ╠═83724a2a-b8ec-4008-a02c-e7ff535f6f87
# ╠═2189aa54-5ec3-4da7-abf1-86633bcc3b3d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
