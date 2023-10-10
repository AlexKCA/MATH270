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
##### Assignment by: 
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

# ╔═╡ fcd2332d-1a3d-4b69-9510-90f57fdf96a3
A = no_zero_pivots(3)

# ╔═╡ 132f5ae7-5cec-41fb-9aff-464c0464e726
md"""
##### 2. Find the ``LU`` factorization of this ``3\times 3`` matrix:
"""

# ╔═╡ 4280220d-8457-4c1a-b458-daf997eedc95
B = zero_pivots(3)

# ╔═╡ b0fc4703-dd42-4e29-8ad5-9bc7a9dcdbbd
md"""
##### 3. Find the ``LU`` factorization of this ``4\times 4`` matrix:
"""

# ╔═╡ 6adab8c9-19c4-4946-96fb-e7642ce398bf
C = no_zero_pivots(4)

# ╔═╡ c77a6304-f389-46b0-8fc4-ea96b815cd1c
md"""
##### 4. Find the ``LU`` factorization of this ``4\times 4`` matrix:
"""

# ╔═╡ c7f64759-f747-4c99-a6bc-e2ffa29d63fc
D = zero_pivots(4)

# ╔═╡ 5080dcf5-992e-4d3b-8409-c9d0f32289cc
md"""
##### **Extra Credit**: Find the ``LDU`` factorization of this ``4\times 4`` matrix, where:
* ``L`` is a lower triangular matrix with ``1``'s on the main diagonal.
* ``D`` is a *diagonal matrix* with numbers on the main diagonal and ``0``'s elsewhere.
* ``U`` is an upper triangular matrix with ``1``'s on the main diagonal.
"""

# ╔═╡ be0404c0-54ce-47f2-b464-586822e970ff
M = no_zero_pivots(4)

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

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "19cc70172396ceef58184fe4a17f4580912aeedf"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

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
version = "5.7.0+0"
"""

# ╔═╡ Cell order:
# ╟─6e51732b-4f10-434f-abd1-1af6ca5e3e11
# ╟─a35ff42e-db2d-47e6-8457-ffa7636677ca
# ╟─bcc7a290-6817-4746-95a0-7ca22a69c424
# ╠═7a54a485-b394-4b18-855f-1b544546f438
# ╟─a187d952-63a9-11ee-39df-8d663e851841
# ╟─78f18405-2cd0-4cbe-844c-895ffd306c18
# ╠═fcd2332d-1a3d-4b69-9510-90f57fdf96a3
# ╟─132f5ae7-5cec-41fb-9aff-464c0464e726
# ╠═4280220d-8457-4c1a-b458-daf997eedc95
# ╟─b0fc4703-dd42-4e29-8ad5-9bc7a9dcdbbd
# ╠═6adab8c9-19c4-4946-96fb-e7642ce398bf
# ╟─c77a6304-f389-46b0-8fc4-ea96b815cd1c
# ╠═c7f64759-f747-4c99-a6bc-e2ffa29d63fc
# ╟─5080dcf5-992e-4d3b-8409-c9d0f32289cc
# ╠═be0404c0-54ce-47f2-b464-586822e970ff
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
