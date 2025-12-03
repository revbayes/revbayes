## name
mvCharacterHistory
## title
Character history move
## description
A proposal to change the (discrete) character history, at a node or along a branch, for a data augmentation-based CTMC distribution (`dnPhyloCTMCDASiteIID` or `dnPhyloCTMCDASequence`).

## details
The `graph` argument specifices where the move is put, the `proposal` argument specifies the type of sampler, while the `type` argument specifies the type of character.

For state-dependent trait evolution (e.g, `dnPhyloOUSD`), moves can be put at the root node (`graph="root"`), at all the internal nodes (`graph="node"`), along branches (`graph="path"`), and at the tips (`graph="tip"`).
In all these cases, a rejection sampler is used (`proposal="rejection"`).
Note that if a move is put at a node, the sampler will propose new character history for the node and its ancestral and descedant branches (if exist), without changing the state of the ancestral and descedant nodes.
If a move is put along a branch, the sampler will propose new character history for the branch, without changing the state of the ancestral and descedant nodes.
The root rejection sampler should only be used when there is only one character.
Furthermore, if some tips have unobserved character states, or if a hidden-state model is used for the CTMC distribution, the character type should be set to natural numbers (`type="NaturalNumbers"`).
Otherwise, the character type should be set to standard (`type="Standard"`).

Refer to respective tutorials for other usage of the data augmentation-based CTMC distribution and the character history move, in particular, biogeographic reconstruction and modeling host repertoire evolution.

## authors
Priscilla Lau
## see_also
## example
The following example shows the set-up of character history moves when the CTMC distribution is used together with a trait evolution model (e.g, `dnPhyloOUSD`).

X ~ dnPhyloCTMCDASiteIID(tree, Q, branchRates=1, type="Standard", nSites=1)

moves.append( mvCharacterHistory(ctmc=X, qmap_site=Q, graph="root",   proposal="rejection", weight=10.0) )
moves.append( mvCharacterHistory(ctmc=X, qmap_site=Q, graph="node",   proposal="rejection", weight=10.0) )
moves.append( mvCharacterHistory(ctmc=X, qmap_site=Q, graph="branch", proposal="rejection", weight=10.0) )
moves.append( mvCharacterHistory(ctmc=X, qmap_site=Q, graph="tip",    proposal="rejection", weight=10.0) )

## references
