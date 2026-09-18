# Instructions

The current repository holds a diversity index called multiplicity as
well as several other indices.

The distance based multiplicity (@R/multiplicity.distance.R) has a
specific implementation called: `multiplicity.distance.by_blocks`. This
was designed for when we have several `block` matrices that correspond
to the distances of subunits inside a unit and we assume that any
distance between subunits of different units is equal to `sigma` (the
value when two subunits or units are considred to be complety distant).
Thus this implementation allows to compute the distance based
multiplicity of an otherwise very large matrix, but where most of the
values are efectively sigma, by only needing the value of the inner unit
subunits (blocks).

I need you to create this type of implementation for the index:
redundancy (@R/other_indices.R)

The scheme is very similar to what we do in
multiplicity.distance.by_blocks. Furthermore,
`multiplicity.distance.by_blocks` computes the necesary Rao’s Q that is
needed in the implementation.

Please review the code, make sure you understand how the block formula
works and give me a plan. Please document the formula as the rest is
documented so that we can export the documentation automatically.
Include also test for this function, were, for small cases the blocks
and non blocks implementation should match.
