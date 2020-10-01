# Create Merrill's Algorithm Function

# The purpose of this script is to create a function that completes
# the recursive calculations of Merrill's algorithm

# ps2, iv2, jout2, labels2 = mCalc( ps1, value1, jout1, iv1, labels1, Nvectors, st, verbose, f)
# mCalc takes in the following values and returns a tuple

# ps1: primitive simplex
# value1: values to label on - in our cases this is ex. demands
# jout1: the column we are going to label
# iv1: the initial vertex
# labels1: labels for primitive simplex
# Nvectors: how many vertices are we searching over
# st: 0/1 is this the first run (if so mCalc initializes the first run)
# verbose: true/false should mCalc write out  to the diagnostic file?
# f: file object for diagnostic file

# ps2, iv2, jout2, labels2 = a tuple representing all these objects
#								as output

# Example:
# v = 2 # number of vertices
# myiv = ones(v)*10 # initial vertex of length v

# if (iter == 1)
#	ps4, iv4, jout4, labels4 = mCalc( 0, zeros(v), 1, myiv, 0, v, 1, false, dfile) # starts primitive simplex
# else
#	ps4, iv4, jout4, labels4 = mCalc( ps4, exDem, jout4, iv4, labels4, v, 0, false, dfile) # continues
# end

function dm( label, value, f)

	write(f, "$(label)\n")

	for i = 1:size(value)[1]
		write(f, "$(value[i,:])")
		write(f, "\n")
	end

end

# Open diagnostic file


function mCalc( ps1, value, jout, iv, labels, Nvectors, st, verbose, f)


	# ------------------------------------------------------------------
	#                  Check Simplex Labels
	# ------------------------------------------------------------------

	# Labeling for vertex on real simplex
	if (ps1[ 1, jout] == 0.) && (st == 0)

		c = 1
		j = 1
		while c == 1

			if(j > Nvectors)
				c = 2
			else
				if value[j] > 0.
					c = 2
					labels[jout] = j
				else
					j = j + 1
				end
			end
		end
	# Print out diagnostic information
	if(verbose)
		write(f, "S1 -> jout: $(jout)\n")
		dm("S1 -> ps:", ps1, f)
		write(f, "S1 -> labels: $(labels)\n")
		write(f, "S1 -> values: $(value)\n")
	end

	end


	# Check for major step
	major = 0
	a = falses(Nvectors+1)
	if st == 1
		major = 1
	else
		# Do we have a full real simplex?
		for i = 1:(Nvectors+1)
			if ps1[1, i] == 0
				a[i] = true
			end
		end
		b = (sum(a) == Nvectors)
		#println("b: ", b)
		#println("a: ", a)
		#println("sum a: ", sum(a))

		if b
			# Have a full real simplex, but
			# is it completely labeled?
			d = size(unique(labels[a]))[1] == Nvectors
		end

		if b && d
		# Set the major indicator
			major = 1

			for i = 1:Nvectors
				iv[i] = ps1[ i+1, jout] * 2  # Refinement Factor
			end

		end

	end

	# ------------------------------------------------------------------
	#                  Major Iteration
	# ------------------------------------------------------------------

	ps3 = zeros( Nvectors+1, Nvectors+1)
	labelsmaj = zeros(Nvectors+1)
	if major == 1

		println("***********************************\n")
		println(" From Merrill -> Reached Major Axis\n")
		println("***********************************\n")


		k = 1.0
		for r = 1:(Nvectors+1)

			if(r > 1)
				println("vertex = $(iv[r-1]) ; value = $(value[r-1])\n")
			end

			for c = 1:(Nvectors+1)

				if c == 1 && r > 1
					ps3[ r, c] = iv[r-1]
				end
				if r == 1
					ps3[ r, c] = 1.
				end
				if (r>1) && (c>1)
					ps3[ r, c] = iv[r-1]
				end
				if (r>1) && (c>1) && (c==r)
					ps3[ r, c] = iv[r-1] - 1
				end

				labelsmaj[c] = c - 1

			end

		end

		# Create Labels
		labelsmaj[1] = 1

		ps3[ 1, 1] = 0

		ps1 = ps3

		jout = 1
		labels = labelsmaj

	end

	# ------------------------------------------------------------------
	#                  Minor Iteration
	# ------------------------------------------------------------------

	# Replacement operation
	x = 1
	t = 1
	myz = Nvectors
	while (x == 1)

		# -------------------------- Artificial Axis Labeling
		c = 2
		j = 2
		# Are we are on the artificial simplex
		if ps1[ 1, jout] == 1.
			c = 1
			#println("Artificial")
		end
		# If we are on the artificial simplex - labeling rule
		while c == 1
			if ps1[ j, jout] < iv[j - 1]
				labels[jout] = j - 1
				c = 2
			else
				if(j > Nvectors)
					c = 2
				end
				j = j + 1
			end
		end

		if (ps1[ 1, jout] == 0) && (myz != Nvectors)

			c = 1
			j = 2
			while c == 1
			# This first loop checks if any of the vertices are at
			# the edge of the simplex - i.e. a coordinate is zero
				if ps1[ j, jout] == 0
					c = 2
					labels[jout] = j - 1
				else
					if(j > Nvectors)
						c = 2
					end
					j = j + 1
				end
			end
		end

		if verbose
			write(f, "S2 -> jout: $(jout)\n")
			dm("S2 -> ps:", ps1, f)
			write(f, "S2 -> labels: $(labels)\n")
		end

		# --------------------------- Find which vertex needs to be replaced
		c = 1
		j = 1
		while(c == 1)
			if (labels[j] == labels[jout]) && (j != jout)
				jout = j
				c = 2
			else
				if(j > Nvectors)
					c = 2
				end
				j = j + 1
			end
		end

		# Now perform the replacement
		j1 = jout - 1
		if(j1 < 1)
			j1 = (Nvectors + 1)
		end
		j2 = jout + 1
		if j2 > (Nvectors + 1)
			j2 = 1
		end

		if verbose
			write(f, "Replacing column $(jout)\n")
		end
		# Replacement operation
		for i = 1:(Nvectors+1)
			ps1[ i, jout] = ps1[ i, j2] + ps1[ i, j1] - ps1[ i, jout]
		end

		myz = sum(ps1[2:(Nvectors+1), jout] .> 0.)
		if verbose
			write(f, "myz = $(myz)\n")
			write(f, "ps1 = $(ps1[1, jout])")
		end

		# Is our new vertex on the artificial simplex?
		if (ps1[ 1, jout] == 0) && (myz == Nvectors)
			x = 2
			labels[jout] = 0
		end

		t = t + 1
		if t > 20
			x = 2
		end

	end

	return ps1, iv, jout, labels

end
