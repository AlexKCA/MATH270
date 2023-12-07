### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ ead6e315-8f6d-45f2-a638-717e40964c43
using Latexify

# ╔═╡ 8440c79a-5999-11ee-3060-ab4b8dd8cd05
md"""
# Elementary Row Operations by Matrices
"""

# ╔═╡ 5e99a294-1fb4-4daa-84ea-6227d68fa331
md"""
### Floating Point and Rational Numbers in Julia

In Julia there are a variety of number types. Among them are Rational number types and Floating point number types.

Floating point number types are decimal numbers with limited precision, whereas Rational number types represent exact rational values.

Floating point numbers (as fractions) are represented with a single `/` division sign, whereas Rational numbers are represented with a `//` division sign.
"""

# ╔═╡ 6c0a20bf-4996-4110-8b95-802b6a33678e
md"""
##### Examples of calculations with Floating point and Rational numbers
"""

# ╔═╡ a5bf79d3-d488-48d2-b988-62f33b6a337e
md"""
* 2/3 = $(2/3)
* 2//3 = $(2//3)
* 2/7 + 3/5 - 4/3 = $(2/7 + 3/5 - 4/3)
* 2//7 + 3//5 - 4//3 = $(2//7 + 3//5 - 4//3)
* (-3/7 + 4/5) / (5/8 - 5/6) = $((-3/7 + 4/5) / (5/8 - 5/6))
* (-3//7 + 4//5) / (5//8 - 5//6) = $((-3//7 + 4//5) / (5//8 - 5//6))
"""

# ╔═╡ 1d96ffc1-e433-4ac6-a56f-8113684fe2a3
md"""
### Matrices in Julia

A matrix in Julia is delimited by brackets [, ] with column entries separated by spaces and rows separated by semicolons ;.

For example, the ``(3\times 4)`` matrix:

$A = \begin{bmatrix} 1 & 2 & -3 & -2\\ -1 & 0 & 2 & 1\\ 3 & -1 & 1 & 3\end{bmatrix}$

is entered in Julia as:
```julia
A = [1 2 -3 -2; -1 0 2 1; 3 -1 1 3]
```
This assigns the matrix to the variable "A".

Since multiple spaces and new lines are ignored in Julia, we can put each row on a new line and add spacing as necessary to align columns to make reading the matrix easier:
```julia
A = 
[1  2 -3 -2;
-1  0  2  1;
 3 -1  1  3]
```
The numbers in this matrix are treated by Julia as Floating point numbers. If we wanted Julia to treat them as Rational numbers, we would enter the matrix as:
```julia
A =
[1//1  2//1 -3//1 -2//1;
-1//1  0//1  2//1  1//1;
 3//1 -1//1  1//1  3//1]
```
where we again used spacing to make the matrix entries easier to read.
"""

# ╔═╡ 65ef861e-1763-4bd3-8597-c1bd5869a171
md"""
# Exercises

In the exercises below, use matrix multiplication to get each matrix into upper triangular form, following my example.
"""

# ╔═╡ d63b665b-d28a-4497-9e69-727f6ba87acb
md"""
### My Example
"""

# ╔═╡ a0656cb7-01db-4f0a-871d-96bd84215958
md"""
First we create a ``(3\times 3)`` matrix of random rational numbers. If the upper left number (the ``(1,1)`` entry) in the matrix is zero, instantiate this cell until that entry is not zero. We give this a "``0``" label.
"""

# ╔═╡ ef590265-d652-44bb-ac8e-817d7e63d6f3
md"""
Now we make a copy of the previous matrix and give it a "``1``" label. After making the copy, disable the cell containing the randomly generated matrix. Disabling that cell will keep the values in that matrix from changing.
"""

# ╔═╡ b8c7ff4a-95d5-452b-9769-572e55e2dea2
A1 = [-4//1 -7//1 -6//1; 6//1 -3//1 4//1; 2//1 -6//1 -1//1]

# ╔═╡ 3f3f7aa9-55ed-450d-89b5-1b998a14d3fe
md"""
We are now ready to start transforming the matrix into upper triangular form. We do this by multiplying the appropriate ``(3\times 3)`` matrices on the left to effectuate row operations. Label each successive transformed matrix with the next succeding number. For example, the transformed ``A1`` matrix becomes matrix ``A2``.
"""

# ╔═╡ edea8fb0-f2aa-4a99-888e-c615b9dc0788
md"""
Here I get a ``1`` in the ``(1,1)`` entry by multiplying the row by ``-1/4``. This makes it easier to eliminate the numbers below it.
"""

# ╔═╡ 6aa85292-91bc-4bfb-a46a-aa909476730c
A2 = [-1//4 0 0; 
	      0 1 0; 
	      0 0 1] * A1

# ╔═╡ a495a177-ff12-435e-a1f3-3111ea651c33
md"""
Next, I use the ``1`` in the ``(1,1)`` position to eliminate the numbers below it. I do this by replacing the second row with ``-6`` times the first row plus the second row, and replacing the third row with ``-2`` times the first row plus the third row.
"""

# ╔═╡ bc225350-ead1-4b92-9752-39b450f05d15
A3 = [1 0 0;
	 -6 1 0;
	 -2 0 1] * A2

# ╔═╡ 7f29f366-fa5f-4a04-8dea-d76e96c3c4bb
md"""
Next, I get a ``1`` in the ``(2,2)`` position by multiplying that row by ``-2/27``. This will make it easier to eliminate the number below it.
"""

# ╔═╡ d2d10361-3387-4f51-a0f8-135ec23c0b46
A4 = [1 0 0;
	  0 -2//27 0;
	  0 0 1] * A3

# ╔═╡ e498cb31-cf72-4ed4-b90e-a1b40e36fc6b
md"""
Finally, I replace the third row with ``19/2`` times the second row plus the third row. This gets the matrix in upper triangular form.
"""

# ╔═╡ ae7f16f4-7378-469e-8604-9fbb82f8ae5e
A5 = [1 0 0;
	  0 1 0;
	  0 19//2 1] * A4

# ╔═╡ a623aadd-0207-4b37-8e81-4e2c2c1a34a9
md"""
Lastly, I express the answer in a more readable form using latexify
"""

# ╔═╡ 5aa5dd9d-f539-4f10-8900-271166c688da
latexify(A5)

# ╔═╡ d5650e88-d0ae-4fbd-a3ae-64a9eacfbf31
md"""
### Your turn
"""

# ╔═╡ 283808e9-a74d-485d-893c-606bbf054a1a
md"""
###### Reduce a ``(3\times 3)`` matrix
"""

# ╔═╡ bd0a7377-3102-463b-81a7-69c11d0ac187
B0 = convert(Matrix{Rational}, rand(-9:9, (3,3)))

# ╔═╡ b3fe58a8-6b46-466b-a099-c15895315f22
md"""
###### Reduce a ``(4\times 4)`` matrix

After you finish reducing the ``(3\times 3)`` matrix to upper triangular form, reduce this ``(4\times 4)`` matrix using the same method as described above.
"""

# ╔═╡ dea1af26-fcb1-48fe-b65d-883e1543566f
C0 = convert(Matrix{Rational}, rand(-9:9, (4,4)))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"

[compat]
Latexify = "~0.16.1"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.0"
manifest_format = "2.0"
project_hash = "4f8d13a58b0f5bd3a860c76bb77c72b934f47f7c"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

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

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9ee1618cbf5240e6d4e0371d6f24065083f60c48"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.11"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

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
"""

# ╔═╡ Cell order:
# ╟─ead6e315-8f6d-45f2-a638-717e40964c43
# ╟─8440c79a-5999-11ee-3060-ab4b8dd8cd05
# ╟─5e99a294-1fb4-4daa-84ea-6227d68fa331
# ╟─6c0a20bf-4996-4110-8b95-802b6a33678e
# ╟─a5bf79d3-d488-48d2-b988-62f33b6a337e
# ╟─1d96ffc1-e433-4ac6-a56f-8113684fe2a3
# ╟─65ef861e-1763-4bd3-8597-c1bd5869a171
# ╟─d63b665b-d28a-4497-9e69-727f6ba87acb
# ╟─a0656cb7-01db-4f0a-871d-96bd84215958
# ╠═72eb5e23-1315-4730-b431-2d2b1b01b14b
# ╟─ef590265-d652-44bb-ac8e-817d7e63d6f3
# ╠═b8c7ff4a-95d5-452b-9769-572e55e2dea2
# ╟─3f3f7aa9-55ed-450d-89b5-1b998a14d3fe
# ╟─edea8fb0-f2aa-4a99-888e-c615b9dc0788
# ╠═6aa85292-91bc-4bfb-a46a-aa909476730c
# ╟─a495a177-ff12-435e-a1f3-3111ea651c33
# ╠═bc225350-ead1-4b92-9752-39b450f05d15
# ╟─7f29f366-fa5f-4a04-8dea-d76e96c3c4bb
# ╠═d2d10361-3387-4f51-a0f8-135ec23c0b46
# ╟─e498cb31-cf72-4ed4-b90e-a1b40e36fc6b
# ╠═ae7f16f4-7378-469e-8604-9fbb82f8ae5e
# ╟─a623aadd-0207-4b37-8e81-4e2c2c1a34a9
# ╠═5aa5dd9d-f539-4f10-8900-271166c688da
# ╟─d5650e88-d0ae-4fbd-a3ae-64a9eacfbf31
# ╟─283808e9-a74d-485d-893c-606bbf054a1a
# ╠═bd0a7377-3102-463b-81a7-69c11d0ac187
# ╟─b3fe58a8-6b46-466b-a099-c15895315f22
# ╠═dea1af26-fcb1-48fe-b65d-883e1543566f
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
