### A Pluto.jl notebook ###
# v0.19.32

using Markdown
using InteractiveUtils

# ╔═╡ a29226c4-3a2b-4ae1-9fb7-f6f75e04c398
using LinearAlgebra, Latexify

# ╔═╡ 7350caea-8d9a-11ee-128f-fdc72c70c829
md"""
# The Church Sequence

The Church sequence ``\{C_n\}`` is the recursive sequence defined by:
```math
C_1 = 0,\; C_2 = 1,\; C_n = C_{n-2} + 2C_{n-1},\quad\text{for } n\geq 3
```
The first few terms of the sequence are:
```math
0, 1, 2, 5, 12, 29, 70, 169, \ldots
```
"""

# ╔═╡ e6debf46-437b-4da4-a7e4-7900c587ff26
md"""
##### Generating ``\{C_n\}`` recursively
"""

# ╔═╡ f0ab11e3-34fc-4bdc-895e-6835c1f5d42b
function chu(n)
	if n == 1
		return 0
	elseif n == 2
		return 1
	else
		chu(n-2) + 2*chu(n-1)
	end
end

# ╔═╡ c6601323-cfdd-404b-abfd-56de3052840c
function chu2(n)
	lst = [0, 1]
	if n == 1
		return [0]
	elseif n == 2
		return [0, 1]
	else
		for i in 1:(n-2)
			append!(lst, 2*lst[end] + lst[end-1])
		end
		return lst
	end
end

# ╔═╡ c9caed98-f184-47d5-840c-75ddfe2be2d5
chu2(15)

# ╔═╡ 7a0f6f76-c816-4353-8318-5d542c55b654
[chu(n) for n in 1:15]

# ╔═╡ 5365e439-cd8c-443d-9793-f154e6f18425
md"""
##### Generating ``\{C_n\}`` using matrix exponentials

Find a matrix ``A`` so that:

```math
A\begin{bmatrix}
	a\\b
\end{bmatrix} = 

\begin{bmatrix}
	b\\a+2b
\end{bmatrix}
```
"""

# ╔═╡ 165b1111-8ca4-42f3-8d3f-872d356d6c4a
# Enter the matrix here using the LaTeX \begin{bmatrix}...\end{bmatrix} environment as shown in the cell above

md"""
!!! answer
```math
A = 
```
"""

# ╔═╡ d3d9092f-e359-4694-a206-54be22de8acf
# Define the matrix A in Julia. Uncomment the line below after you enter the matrix, then run the cell.

# A = 

# ╔═╡ db7cbe5c-c91d-43a7-bf9f-4fd73d0e4c58
md"""
Which Church numbers are in the vector ``A^3\begin{bmatrix}0\\1\end{bmatrix}``?
"""

# ╔═╡ 649d54a9-9fa4-46d8-8f67-4442010db3aa
# uncomment the line below, then run the cell. Enter in your answers in the cell below.

# A^3 * [0; 1]

# ╔═╡ 5db735fe-57ae-495c-adb1-1ca89fe1988a
md"""
!!! answer
The Church numbers in the vector ``A^3\begin{bmatrix} 0\\1 \end{bmatrix}`` are:

"""

# ╔═╡ 1fda3b8d-46da-4254-a2ec-671665d70f64
md"""
What does the exponent ``m`` on the matrix ``A`` have to be, so that 
```math
A^m \begin{bmatrix} 0\\1 \end{bmatrix} = \begin{bmatrix} C_{n-1}\\C_n \end{bmatrix}\\?
```
"""

# ╔═╡ b6efe442-1db1-4849-aa2b-d380c6ca3443
# Put the answer here.
md"""
!!! answer
``m =  ``
"""

# ╔═╡ e73e0c9f-cc1c-4047-aefc-ffdc9dbc7a50
md"""
!!! answer
Use a list comprehension using powers of the matrix ``A`` to generate the third to the fifteenth Church numbers
"""

# ╔═╡ f175506e-cb90-4adb-81ac-292ac16ae9d3
# Put the list comprehension here. Verify that the numbers agree with the third to fifteenth Church numbers.



# ╔═╡ 678957e1-509c-4352-8dfd-9ed4ff4dce0c
md"""
#### Diagonalizing the matrix ``A``
"""

# ╔═╡ f62f064f-882b-4638-ab9f-d93327b291f8
md"""
###### Find the eigenvalues and corresponding eigenvectors of the matrix ``A``
"""

# ╔═╡ 6f563de2-ed6a-4a94-b2c8-fb6bfd86f4aa
md"""
!!! answer
Enter the eigenvalues in ``\LaTeX`` below:

```math
\begin{align}
	\lambda_1 &=    \\
	\lambda_2 &=  
\end{align}
```
"""

# ╔═╡ de17315a-8d8e-4b4d-b857-1a553f18c0ab
md"""
!!! answer
Enter the corresponding eigenvectors in ``\LaTeX`` below:

```math
\begin{align}
	v_1 &= \begin{bmatrix}  \end{bmatrix}\\[1.5ex]
	v_2 &= \begin{bmatrix}  \end{bmatrix}
\end{align}
```
"""

# ╔═╡ 06a4d426-ea92-45d8-8b9e-3782dee98d21
md"""
!!! answer
Enter the matrix ``X`` of eigenvectors in ``\LaTeX`` below:

```math
X = \begin{bmatrix}

\end{bmatrix}
```
"""

# ╔═╡ 9b94cb30-5d35-4889-bf6f-e9a729044cac
md"""
!!! answer
Enter the inverse of matrix ``X`` in ``\LaTeX`` below:

```math
X^{-1} = \begin{bmatrix}

\end{bmatrix}
```
"""

# ╔═╡ b10cd556-3866-448a-bbcd-8b20f45217d8
md"""
!!! answer
Enter the diagonal matrix ``\Lambda`` of eigenvalues of matrix ``A`` in ``\LaTeX`` below:

```math
\Lambda = \begin{bmatrix}

\end{bmatrix}
```
"""

# ╔═╡ bc9c3361-2dcf-4d3c-9e5b-abaec0eaa28b
md"""
!!! answer
Using the diagonalization ``A = X\Lambda X^{-1}``, find entries in the vector ``A^m \begin{bmatrix} 0\\1 \end{bmatrix}``. Enter your answer below in ``\LaTeX``.

```math
A^m \begin{bmatrix} 0\\1 \end{bmatrix} = 
\begin{bmatrix}

\end{bmatrix}
```
"""

# ╔═╡ f6744705-e888-4198-81fa-849aa739febc
md"""
!!! answer
Based on your previous work, what is a formula for the ``n``th Church number? Enter your answer in ``\LaTeX`` below:

```math
C_n = 
```
"""

# ╔═╡ f8c33c57-4334-4592-9b27-c84d8227367d
md"""
!!! answer
Using the formula you discovered for ``C_n``, make a julia function C(n) for the ``n``th Church number. Enter it in the cell below.
"""

# ╔═╡ 2eef3fd6-d386-495c-96dc-50dde6e072a9
# uncomment the line below and run the cell after you enter in your formula for C(n)

# C(n) = 

# ╔═╡ de1b9a28-b38d-46eb-9019-11f967e72841
md"""
!!! answer
Verify that your function C(n) produces the correct results by making a list comprehension for C(3) to C(15) and checking that the results agree with the third to fifteenth Church numbers. Enter the list comprehension in the cell below.
"""

# ╔═╡ 5c5d64d2-47cf-4c17-9c1e-45932116ce37
# Put your list comprehension here, then run the cell



# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[compat]
Latexify = "~0.16.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "3ed0b82707adb6cd636bc44fc0357f595e3c8af1"

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
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

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
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

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
# ╟─a29226c4-3a2b-4ae1-9fb7-f6f75e04c398
# ╟─7350caea-8d9a-11ee-128f-fdc72c70c829
# ╟─e6debf46-437b-4da4-a7e4-7900c587ff26
# ╠═f0ab11e3-34fc-4bdc-895e-6835c1f5d42b
# ╠═c6601323-cfdd-404b-abfd-56de3052840c
# ╠═c9caed98-f184-47d5-840c-75ddfe2be2d5
# ╠═7a0f6f76-c816-4353-8318-5d542c55b654
# ╠═5365e439-cd8c-443d-9793-f154e6f18425
# ╠═165b1111-8ca4-42f3-8d3f-872d356d6c4a
# ╠═d3d9092f-e359-4694-a206-54be22de8acf
# ╟─db7cbe5c-c91d-43a7-bf9f-4fd73d0e4c58
# ╠═649d54a9-9fa4-46d8-8f67-4442010db3aa
# ╠═5db735fe-57ae-495c-adb1-1ca89fe1988a
# ╟─1fda3b8d-46da-4254-a2ec-671665d70f64
# ╠═b6efe442-1db1-4849-aa2b-d380c6ca3443
# ╟─e73e0c9f-cc1c-4047-aefc-ffdc9dbc7a50
# ╠═f175506e-cb90-4adb-81ac-292ac16ae9d3
# ╟─678957e1-509c-4352-8dfd-9ed4ff4dce0c
# ╟─f62f064f-882b-4638-ab9f-d93327b291f8
# ╠═6f563de2-ed6a-4a94-b2c8-fb6bfd86f4aa
# ╠═de17315a-8d8e-4b4d-b857-1a553f18c0ab
# ╠═06a4d426-ea92-45d8-8b9e-3782dee98d21
# ╠═9b94cb30-5d35-4889-bf6f-e9a729044cac
# ╠═b10cd556-3866-448a-bbcd-8b20f45217d8
# ╠═bc9c3361-2dcf-4d3c-9e5b-abaec0eaa28b
# ╠═f6744705-e888-4198-81fa-849aa739febc
# ╟─f8c33c57-4334-4592-9b27-c84d8227367d
# ╠═2eef3fd6-d386-495c-96dc-50dde6e072a9
# ╟─de1b9a28-b38d-46eb-9019-11f967e72841
# ╠═5c5d64d2-47cf-4c17-9c1e-45932116ce37
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
