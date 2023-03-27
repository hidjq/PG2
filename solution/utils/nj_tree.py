import numpy as np
from six import StringIO
from six import StringIO, string_types
from skbio.tree import TreeNode
class DissimilarityMatrix(object):
    """Store dissimilarities between objects.
    A `DissimilarityMatrix` instance stores a square, hollow, two-dimensional
    matrix of dissimilarities between objects. Objects could be, for example,
    samples or DNA sequences. A sequence of IDs accompanies the
    dissimilarities.
    Methods are provided to load and save dissimilarity matrices from/to disk,
    as well as perform common operations such as extracting dissimilarities
    based on object ID.
    Parameters
    ----------
    data : array_like or DissimilarityMatrix
        Square, hollow, two-dimensional ``numpy.ndarray`` of dissimilarities
        (floats), or a structure that can be converted to a ``numpy.ndarray``
        using ``numpy.asarray``. Can instead be a `DissimilarityMatrix` (or
        subclass) instance, in which case the instance's data will be used.
        Data will be converted to a float ``dtype`` if necessary. A copy will
        *not* be made if already a ``numpy.ndarray`` with a float ``dtype``.
    ids : sequence of str, optional
        Sequence of strings to be used as object IDs. Must match the number of
        rows/cols in `data`. If ``None`` (the default), IDs will be
        monotonically-increasing integers cast as strings, with numbering
        starting from zero, e.g., ``('0', '1', '2', '3', ...)``.
    Attributes
    ----------
    data
    ids
    dtype
    shape
    size
    T
    See Also
    --------
    DistanceMatrix
    Notes
    -----
    The dissimilarities are stored in redundant (square-form) format [1]_.
    The data are not checked for symmetry, nor guaranteed/assumed to be
    symmetric.
    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
    """
    default_write_format = 'lsmat'
    # Used in __str__
    _matrix_element_name = 'dissimilarity'

    @classmethod
    def from_file(cls, lsmat_f, delimiter='\t'):
        """Load dissimilarity matrix from delimited text file.
        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``from_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``read``, which is a more general method for deserializing
           dissimilarity/distance matrices. ``read`` supports multiple file
           formats, automatic file format detection, etc. by taking advantage
           of scikit-bio's I/O registry system. See :mod:`skbio.io` for more
           details.
        Creates a ``DissimilarityMatrix`` (or subclass) instance from a
        ``lsmat`` formatted file. See :mod:`skbio.io.lsmat` for the format
        specification.
        Parameters
        ----------
        lsmat_f: filepath or filehandle
            File to read from.
        delimiter : str, optional
            String delimiting elements in `lsmat_f`.
        Returns
        -------
        DissimilarityMatrix
            Instance of type `cls` containing the parsed contents of `lsmat_f`.
        See Also
        --------
        read
        """
        warnings.warn(
            "DissimilarityMatrix.from_file and DistanceMatrix.from_file are "
            "deprecated and will be removed in scikit-bio 0.3.0. Please "
            "update your code to use DissimilarityMatrix.read and "
            "DistanceMatrix.read.", UserWarning)
        return cls.read(lsmat_f, format='lsmat', delimiter=delimiter)

    def to_file(self, out_f, delimiter='\t'):
        """Save dissimilarity matrix to file as delimited text.
        .. note:: Deprecated in scikit-bio 0.2.0-dev
           ``to_file`` will be removed in scikit-bio 0.3.0. It is replaced by
           ``write``, which is a more general method for serializing
           dissimilarity/distance matrices. ``write`` supports multiple file
           formats by taking advantage of scikit-bio's I/O registry system.
           See :mod:`skbio.io` for more details.
        Serializes dissimilarity matrix as a ``lsmat`` formatted file. See
        :mod:`skbio.io.lsmat` for the format specification.
        Parameters
        ----------
        out_f : filepath or filehandle
            File to write to.
        delimiter : str, optional
            Delimiter used to separate elements in output format.
        See Also
        --------
        write
        """
        warnings.warn(
            "DissimilarityMatrix.to_file and DistanceMatrix.to_file are "
            "deprecated and will be removed in scikit-bio 0.3.0. Please "
            "update your code to use DissimilarityMatrix.write and "
            "DistanceMatrix.write.", UserWarning)
        self.write(out_f, format='lsmat', delimiter=delimiter)

    def __init__(self, data, ids=None):
        if isinstance(data, DissimilarityMatrix):
            data = data.data
        data = np.asarray(data, dtype='float')

        if ids is None:
            ids = (str(i) for i in range(data.shape[0]))
        ids = tuple(ids)

        self._validate(data, ids)

        self._data = data
        self._ids = ids
        self._id_index = self._index_list(self._ids)

    @property
    def data(self):
        """Array of dissimilarities.
        A square, hollow, two-dimensional ``numpy.ndarray`` of dissimilarities
        (floats). A copy is *not* returned.
        Notes
        -----
        This property is not writeable.
        """
        return self._data

    @property
    def ids(self):
        """Tuple of object IDs.
        A tuple of strings, one for each object in the dissimilarity matrix.
        Notes
        -----
        This property is writeable, but the number of new IDs must match the
        number of objects in `data`.
        """
        return self._ids

    @ids.setter
    def ids(self, ids_):
        ids_ = tuple(ids_)
        self._validate(self.data, ids_)
        self._ids = ids_
        self._id_index = self._index_list(self._ids)

    @property
    def dtype(self):
        """Data type of the dissimilarities."""
        return self.data.dtype

    @property
    def shape(self):
        """Two-element tuple containing the dissimilarity matrix dimensions.
        Notes
        -----
        As the dissimilarity matrix is guaranteed to be square, both tuple
        entries will always be equal.
        """
        return self.data.shape

    @property
    def size(self):
        """Total number of elements in the dissimilarity matrix.
        Notes
        -----
        Equivalent to ``self.shape[0] * self.shape[1]``.
        """
        return self.data.size

    @property
    def T(self):
        """Transpose of the dissimilarity matrix.
        See Also
        --------
        transpose
        """
        return self.transpose()

    def transpose(self):
        """Return the transpose of the dissimilarity matrix.
        Notes
        -----
        A deep copy is returned.
        Returns
        -------
        DissimilarityMatrix
            Transpose of the dissimilarity matrix. Will be the same type as
            `self`.
        """
        return self.__class__(self.data.T.copy(), deepcopy(self.ids))

    def index(self, lookup_id):
        """Return the index of the specified ID.
        Parameters
        ----------
        lookup_id : str
            ID whose index will be returned.
        Returns
        -------
        int
            Row/column index of `lookup_id`.
        Raises
        ------
        MissingIDError
            If `lookup_id` is not in the dissimilarity matrix.
        """
        if lookup_id in self:
            return self._id_index[lookup_id]
        else:
            raise MissingIDError(lookup_id)

    def redundant_form(self):
        """Return an array of dissimilarities in redundant format.
        As this is the native format that the dissimilarities are stored in,
        this is simply an alias for `data`.
        Returns
        -------
        ndarray
            Two-dimensional ``numpy.ndarray`` of dissimilarities in redundant
            format.
        Notes
        -----
        Redundant format is described in [1]_.
        Does *not* return a copy of the data.
        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
        """
        return self.data

    def copy(self):
        """Return a deep copy of the dissimilarity matrix.
        Returns
        -------
        DissimilarityMatrix
            Deep copy of the dissimilarity matrix. Will be the same type as
            `self`.
        """
        # We deepcopy IDs in case the tuple contains mutable objects at some
        # point in the future.
        return self.__class__(self.data.copy(), deepcopy(self.ids))

    def filter(self, ids, strict=True):
        """Filter the dissimilarity matrix by IDs.
        Parameters
        ----------
        ids : iterable of str
            IDs to retain. May not contain duplicates or be empty. Each ID must
            be present in the dissimilarity matrix.
        strict : bool, optional
            If `strict` is ``True`` and an ID that is not found in the distance
            matrix is found in `ids`, a ``MissingIDError`` exception will be
            raised, otherwise the ID will be ignored.
        Returns
        -------
        DissimilarityMatrix
            Filtered dissimilarity matrix containing only the IDs specified in
            `ids`. IDs will be in the same order as they appear in `ids`.
        Raises
        ------
        MissingIDError
            If an ID in `ids` is not in the object's list of IDs.
        """
        if strict:
            idxs = [self.index(id_) for id_ in ids]
        else:
            # get the indices to slice the inner numpy array
            idxs = []
            # save the IDs that were found in the distance matrix
            found_ids = []
            for id_ in ids:
                try:
                    idxs.append(self.index(id_))
                    found_ids.append(id_)
                except MissingIDError:
                    pass
            ids = found_ids

        filtered_data = self._data[idxs][:, idxs]
        return self.__class__(filtered_data, ids)

    def __str__(self):
        """Return a string representation of the dissimilarity matrix.
        Summary includes matrix dimensions, a (truncated) list of IDs, and
        (truncated) array of dissimilarities.
        Returns
        -------
        str
            String representation of the dissimilarity matrix.
        .. shownumpydoc
        """
        return '%dx%d %s matrix\nIDs:\n%s\nData:\n' % (
            self.shape[0], self.shape[1], self._matrix_element_name,
            self._pprint_ids()) + str(self.data)

    def __eq__(self, other):
        """Compare this dissimilarity matrix to another for equality.
        Two dissimilarity matrices are equal if they have the same shape, IDs
        (in the same order!), and have data arrays that are equal.
        Checks are *not* performed to ensure that `other` is a
        `DissimilarityMatrix` instance.
        Parameters
        ----------
        other : DissimilarityMatrix
            Dissimilarity matrix to compare to for equality.
        Returns
        -------
        bool
            ``True`` if `self` is equal to `other`, ``False`` otherwise.
        .. shownumpydoc
        """
        equal = True

        # The order these checks are performed in is important to be as
        # efficient as possible. The check for shape equality is not strictly
        # necessary as it should be taken care of in np.array_equal, but I'd
        # rather explicitly bail before comparing IDs or data. Use array_equal
        # instead of (a == b).all() because of this issue:
        #     http://stackoverflow.com/a/10582030
        try:
            if self.shape != other.shape:
                equal = False
            elif self.ids != other.ids:
                equal = False
            elif not np.array_equal(self.data, other.data):
                equal = False
        except AttributeError:
            equal = False

        return equal

    def __ne__(self, other):
        """Determine whether two dissimilarity matrices are not equal.
        Parameters
        ----------
        other : DissimilarityMatrix
            Dissimilarity matrix to compare to.
        Returns
        -------
        bool
            ``True`` if `self` is not equal to `other`, ``False`` otherwise.
        See Also
        --------
        __eq__
        .. shownumpydoc
        """
        return not self == other

    def __contains__(self, lookup_id):
        """Check if the specified ID is in the dissimilarity matrix.
        Parameters
        ----------
        lookup_id : str
            ID to search for.
        Returns
        -------
        bool
            ``True`` if `lookup_id` is in the dissimilarity matrix, ``False``
            otherwise.
        See Also
        --------
        index
        .. shownumpydoc
        """
        return lookup_id in self._id_index

    def __getitem__(self, index):
        """Slice into dissimilarity data by object ID or numpy indexing.
        Extracts data from the dissimilarity matrix by object ID, a pair of
        IDs, or numpy indexing/slicing.
        Parameters
        ----------
        index : str, two-tuple of str, or numpy index
            `index` can be one of the following forms: an ID, a pair of IDs, or
            a numpy index.
            If `index` is a string, it is assumed to be an ID and a
            ``numpy.ndarray`` row vector is returned for the corresponding ID.
            Note that the ID's row of dissimilarities is returned, *not* its
            column. If the matrix is symmetric, the two will be identical, but
            this makes a difference if the matrix is asymmetric.
            If `index` is a two-tuple of strings, each string is assumed to be
            an ID and the corresponding matrix element is returned that
            represents the dissimilarity between the two IDs. Note that the
            order of lookup by ID pair matters if the matrix is asymmetric: the
            first ID will be used to look up the row, and the second ID will be
            used to look up the column. Thus, ``dm['a', 'b']`` may not be the
            same as ``dm['b', 'a']`` if the matrix is asymmetric.
            Otherwise, `index` will be passed through to
            ``DissimilarityMatrix.data.__getitem__``, allowing for standard
            indexing of a ``numpy.ndarray`` (e.g., slicing).
        Returns
        -------
        ndarray or scalar
            Indexed data, where return type depends on the form of `index` (see
            description of `index` for more details).
        Raises
        ------
        MissingIDError
            If the ID(s) specified in `index` are not in the dissimilarity
            matrix.
        Notes
        -----
        The lookup based on ID(s) is quick.
        .. shownumpydoc
        """
        if isinstance(index, string_types):
            return self.data[self.index(index)]
        elif self._is_id_pair(index):
            return self.data[self.index(index[0]), self.index(index[1])]
        else:
            return self.data.__getitem__(index)

    def _validate(self, data, ids):
        """Validate the data array and IDs.
        Checks that the data is at least 1x1 in size, 2D, square, hollow, and
        contains only floats. Also checks that IDs are unique and that the
        number of IDs matches the number of rows/cols in the data array.
        Subclasses can override this method to perform different/more specific
        validation (e.g., see `DistanceMatrix`).
        Notes
        -----
        Accepts arguments instead of inspecting instance attributes to avoid
        creating an invalid dissimilarity matrix before raising an error.
        Otherwise, the invalid dissimilarity matrix could be used after the
        exception is caught and handled.
        """
        num_ids = len(ids)

        if 0 in data.shape:
            raise DissimilarityMatrixError("Data must be at least 1x1 in "
                                           "size.")
        elif len(data.shape) != 2:
            raise DissimilarityMatrixError("Data must have exactly two "
                                           "dimensions.")
        elif data.shape[0] != data.shape[1]:
            raise DissimilarityMatrixError("Data must be square (i.e., have "
                                           "the same number of rows and "
                                           "columns).")
        elif data.dtype != np.double:
            raise DissimilarityMatrixError("Data must contain only floating "
                                           "point values.")
        elif np.trace(data) != 0:
            raise DissimilarityMatrixError("Data must be hollow (i.e., the "
                                           "diagonal can only contain zeros).")
        elif num_ids != len(set(ids)):
            raise DissimilarityMatrixError("IDs must be unique.")
        elif num_ids != data.shape[0]:
            raise DissimilarityMatrixError("The number of IDs must match the "
                                           "number of rows/columns in the "
                                           "data.")

    def _index_list(self, list_):
        return {id_: idx for idx, id_ in enumerate(list_)}

    def _is_id_pair(self, index):
        return (isinstance(index, tuple) and
                len(index) == 2 and
                all(map(lambda e: isinstance(e, string_types), index)))

    def _pprint_ids(self, max_chars=80, delimiter=', ', suffix='...',):
        # Adapted from http://stackoverflow.com/a/250373
        ids_str = delimiter.join(self.ids)

        if len(ids_str) > max_chars:
            truncated = ids_str[:max_chars + 1].split(delimiter)[0:-1]
            ids_str = delimiter.join(truncated) + delimiter + suffix

        return ids_str


class DissimilarityMatrixError(Exception):
    """General error for dissimilarity matrix validation failures."""
    pass


class DistanceMatrixError(DissimilarityMatrixError):
    """General error for distance matrix validation failures."""
    pass


class MissingIDError(DissimilarityMatrixError):
    """Error for ID lookup that doesn't exist in the dissimilarity matrix."""

    def __init__(self, missing_id):
        super(MissingIDError, self).__init__()
        self.args = ("The ID '%s' is not in the dissimilarity matrix." %
                     missing_id,)


class DistanceMatrix(DissimilarityMatrix):
    """Store distances between objects.
    A `DistanceMatrix` is a `DissimilarityMatrix` with the additional
    requirement that the matrix data is symmetric. There are additional methods
    made available that take advantage of this symmetry.
    See Also
    --------
    DissimilarityMatrix
    Notes
    -----
    The distances are stored in redundant (square-form) format [1]_. To
    facilitate use with other scientific Python routines (e.g., scipy), the
    distances can be retrieved in condensed (vector-form) format using
    `condensed_form`.
    `DistanceMatrix` only requires that the distances it stores are symmetric.
    Checks are *not* performed to ensure the other three metric properties
    hold (non-negativity, identity of indiscernibles, and triangle inequality)
    [2]_. Thus, a `DistanceMatrix` instance can store distances that are not
    metric.
    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
    .. [2] http://planetmath.org/metricspace
    """

    # Override here, used in superclass __str__
    _matrix_element_name = 'distance'

    def condensed_form(self):
        """Return an array of distances in condensed format.
        Returns
        -------
        ndarray
            One-dimensional ``numpy.ndarray`` of distances in condensed format.
        Notes
        -----
        Condensed format is described in [1]_.
        The conversion is not a constant-time operation, though it should be
        relatively quick to perform.
        References
        ----------
        .. [1] http://docs.scipy.org/doc/scipy/reference/spatial.distance.html
        """
        return squareform(self._data, force='tovector', checks=False)

    def permute(self, condensed=False):
        """Randomly permute both rows and columns in the matrix.
        Randomly permutes the ordering of rows and columns in the matrix. The
        same permutation is applied to both rows and columns in order to
        maintain symmetry and hollowness. Only the rows/columns in the distance
        matrix are permuted; the IDs are *not* permuted.
        Parameters
        ----------
        condensed : bool, optional
            If ``True``, return the permuted distance matrix in condensed
            format. Otherwise, return the permuted distance matrix as a new
            ``DistanceMatrix`` instance.
        Returns
        -------
        DistanceMatrix or ndarray
            Permuted distances as a new ``DistanceMatrix`` or as a ``ndarray``
            in condensed format.
        See Also
        --------
        condensed_form
        Notes
        -----
        This method does not modify the distance matrix that it is called on.
        It is more efficient to pass ``condensed=True`` than permuting the
        distance matrix and then converting to condensed format.
        """
        order = np.random.permutation(self.shape[0])
        permuted = self._data[order][:, order]

        if condensed:
            return squareform(permuted, force='tovector', checks=False)
        else:
            return self.__class__(permuted, self.ids)

    def _validate(self, data, ids):
        """Validate the data array and IDs.
        Overrides the superclass `_validate`. Performs a check for symmetry in
        addition to the checks performed in the superclass.
        """
        super(DistanceMatrix, self)._validate(data, ids)

        if (data.T != data).any():
            raise DistanceMatrixError("Data must be symmetric.")


def randdm(num_objects, ids=None, constructor=None, random_fn=None):
    """Generate a distance matrix populated with random distances.
    Using the default `random_fn`, distances are randomly drawn from a uniform
    distribution over ``[0, 1)``.
    Regardless of `random_fn`, the resulting distance matrix is guaranteed to
    be symmetric and hollow.
    Parameters
    ----------
    num_objects : int
        The number of objects in the resulting distance matrix. For example, if
        `num_objects` is 3, a 3x3 distance matrix will be returned.
    ids : sequence of str or None, optional
        A sequence of strings to be used as IDs. ``len(ids)`` must be equal to
        `num_objects`. If not provided, IDs will be monotonically-increasing
        integers cast as strings (numbering starts at 1). For example,
        ``('1', '2', '3')``.
    constructor : type, optional
        `DissimilarityMatrix` or subclass constructor to use when creating the
        random distance matrix. The returned distance matrix will be of this
        type. If ``None`` (the default), a `DistanceMatrix` instance will be
        returned.
    random_fn : function, optional
        Function to generate random values. `random_fn` must accept two
        arguments (number of rows and number of columns) and return a 2D
        ``numpy.ndarray`` of floats (or something that can be cast to float).
        If ``None`` (the default), ``numpy.random.rand`` will be used.
    Returns
    -------
    DissimilarityMatrix
        `DissimilarityMatrix` (or subclass) instance of random distances. Type
        depends on `constructor`.
    See Also
    --------
    numpy.random.rand
    """
    if constructor is None:
        constructor = DistanceMatrix
    if random_fn is None:
        random_fn = np.random.rand

    data = np.tril(random_fn(num_objects, num_objects), -1)
    data += data.T

    if not ids:
        ids = map(str, range(1, num_objects + 1))

    return constructor(data, ids)


def nj(dm, disallow_negative_branch_length=True, result_constructor=None):
    """ Apply neighbor joining for phylogenetic reconstruction.
    Parameters
    ----------
    dm : skbio.DistanceMatrix
        Input distance matrix containing distances between OTUs.
    disallow_negative_branch_length : bool, optional
        Neighbor joining can result in negative branch lengths, which don't
        make sense in an evolutionary context. If `True`, negative branch
        lengths will be returned as zero, a common strategy for handling this
        issue that was proposed by the original developers of the algorithm.
    result_constructor : function, optional
        Function to apply to construct the result object. This must take a
        newick-formatted string as input. The result of applying this function
        to a newick-formatted string will be returned from this function. This
        defaults to ``lambda x: TreeNode.read(StringIO(x), format='newick')``.
    Returns
    -------
    TreeNode
        By default, the result object is a `TreeNode`, though this can be
        overridden by passing `result_constructor`.
    See Also
    --------
    TreeNode.root_at_midpoint
    Notes
    -----
    Neighbor joining was initially described in Saitou and Nei (1987) [1]_. The
    example presented here is derived from the Wikipedia page on neighbor
    joining [2]_. The Phylip manual also describes the method [3]_ and Phylip
    itself provides an implementation which is useful for comparison.
    Neighbor joining, by definition, creates unrooted trees. One strategy for
    rooting the resulting trees is midpoint rooting, which is accessible as
    ``TreeNode.root_at_midpoint``.
    References
    ----------
    .. [1] Saitou N, and Nei M. (1987) "The neighbor-joining method: a new
       method for reconstructing phylogenetic trees." Molecular Biology and
       Evolution. PMID: 3447015.
    .. [2] http://en.wikipedia.org/wiki/Neighbour_joining
    .. [3] http://evolution.genetics.washington.edu/phylip/doc/neighbor.html
    Examples
    --------
    Define a new distance matrix object describing the distances between five
    OTUs: a, b, c, d, and e.
    >>> from skbio import DistanceMatrix
    >>> from skbio.tree import nj
    >>> data = [[0,  5,  9,  9,  8],
    ...         [5,  0, 10, 10,  9],
    ...         [9, 10,  0,  8,  7],
    ...         [9, 10,  8,  0,  3],
    ...         [8,  9,  7,  3,  0]]
    >>> ids = list('abcde')
    >>> dm = DistanceMatrix(data, ids)
    Contstruct the neighbor joining tree representing the relationship between
    those OTUs. This is returned as a TreeNode object.
    >>> tree = nj(dm)
    >>> print(tree.ascii_art())
              /-d
             |
             |          /-c
             |---------|
    ---------|         |          /-b
             |          \--------|
             |                    \-a
             |
              \-e
    Again, construct the neighbor joining tree, but instead return the newick
    string representing the tree, rather than the TreeNode object. (Note that
    in this example the string output is truncated when printed to facilitate
    rendering.)
    >>> newick_str = nj(dm, result_constructor=str)
    >>> print(newick_str[:55], "...")
    (d:2.000000, (c:4.000000, (b:3.000000, a:2.000000):3.00 ...
    """
    if dm.shape[0] < 3:
        raise ValueError(
            "Distance matrix must be at least 3x3 to "
            "generate a neighbor joining tree.")

    if result_constructor is None:
        result_constructor = \
            lambda x: TreeNode.read(StringIO(x), format='newick')

    # initialize variables
    node_definition = None

    # while there are still more than three distances in the distance matrix,
    # join neighboring nodes.
    while(dm.shape[0] > 3):
        # compute the Q matrix
        q = _compute_q(dm)

        # identify the pair of nodes that have the lowest Q value. if multiple
        # pairs have equally low Q values, the first pair identified (closest
        # to the top-left of the matrix) will be chosen. these will be joined
        # in the current node.
        idx1, idx2 = _lowest_index(q)
        pair_member_1 = dm.ids[idx1]
        pair_member_2 = dm.ids[idx2]
        # determine the distance of each node to the new node connecting them.
        pair_member_1_len, pair_member_2_len = _pair_members_to_new_node(
            dm, idx1, idx2, disallow_negative_branch_length)
        # define the new node in newick style
        node_definition = "(%s:%f, %s:%f)" % (pair_member_1,
                                              pair_member_1_len,
                                              pair_member_2,
                                              pair_member_2_len)
        # compute the new distance matrix, which will contain distances of all
        # other nodes to this new node
        dm = _compute_collapsed_dm(
            dm, pair_member_1, pair_member_2,
            disallow_negative_branch_length=disallow_negative_branch_length,
            new_node_id=node_definition)

    # When there are three distances left in the distance matrix, we have a
    # fully defined tree. The last node is internal, and its distances are
    # defined by these last three values.
    # First determine the distance between the last two nodes to be joined in
    # a pair...
    pair_member_1 = dm.ids[1]
    pair_member_2 = dm.ids[2]
    pair_member_1_len, pair_member_2_len = \
        _pair_members_to_new_node(dm, pair_member_1, pair_member_2,
                                  disallow_negative_branch_length)
    # ...then determine their distance to the other remaining node, but first
    # handle the trival case where the input dm was only 3 x 3
    node_definition = node_definition or dm.ids[0]
    internal_len = _otu_to_new_node(
        dm, pair_member_1, pair_member_2, node_definition,
        disallow_negative_branch_length=disallow_negative_branch_length)
    # ...and finally create the newick string describing the whole tree.
    newick = "(%s:%f, %s:%f, %s:%f);" % (pair_member_1, pair_member_1_len,
                                         node_definition, internal_len,
                                         pair_member_2, pair_member_2_len)

    # package the result as requested by the user and return it.
    return result_constructor(newick)

def _compute_q(dm):
    """Compute Q matrix, used to identify the next pair of nodes to join.
    """
    q = np.zeros(dm.shape)
    n = dm.shape[0]
    for i in range(n):
        for j in range(i):
            q[i, j] = q[j, i] = \
                ((n - 2) * dm[i, j]) - dm[i].sum() - dm[j].sum()
    return DistanceMatrix(q, dm.ids)
def _lowest_index(dm):
    """Return the index of the lowest value in the input distance matrix.
    If there are ties for the lowest value, the index of top-left most
    occurrence of that value will be returned.
    This should be ultimately be replaced with a new DistanceMatrix object
    method (#228).
    """
    lowest_value = np.inf
    for i in range(dm.shape[0]):
        for j in range(i):
            curr_index = i, j
            curr_value = dm[curr_index]
            if curr_value < lowest_value:
                lowest_value = curr_value
                result = curr_index
    return result
def _compute_collapsed_dm(dm, i, j, disallow_negative_branch_length,
                          new_node_id):
    """Return the distance matrix resulting from joining ids i and j in a node.
    If the input distance matrix has shape ``(n, n)``, the result will have
    shape ``(n-1, n-1)`` as the ids `i` and `j` are collapsed to a single new
    ids.
    """
    in_n = dm.shape[0]
    out_n = in_n - 1
    out_ids = [new_node_id]
    out_ids.extend([e for e in dm.ids if e not in (i, j)])
    result = np.zeros((out_n, out_n))
    for idx1, out_id1 in enumerate(out_ids[1:]):
        result[0, idx1 + 1] = result[idx1 + 1, 0] = _otu_to_new_node(
            dm, i, j, out_id1, disallow_negative_branch_length)
        for idx2, out_id2 in enumerate(out_ids[1:idx1+1]):
            result[idx1+1, idx2+1] = result[idx2+1, idx1+1] = \
                dm[out_id1, out_id2]
    return DistanceMatrix(result, out_ids)


def _otu_to_new_node(dm, i, j, k, disallow_negative_branch_length):
    """Return the distance between a new node and some other node.
    Parameters
    ----------
    dm : skbio.DistanceMatrix
        The input distance matrix.
    i, j : str
        Identifiers of entries in the distance matrix to be collapsed. These
        get collapsed to a new node, internally represented as `u`.
    k : str
        Identifier of the entry in the distance matrix for which distance to
        `u` will be computed.
    disallow_negative_branch_length : bool
        Neighbor joining can result in negative branch lengths, which don't
        make sense in an evolutionary context. If `True`, negative branch
        lengths will be returned as zero, a common strategy for handling this
        issue that was proposed by the original developers of the algorithm.
    """
    k_to_u = 0.5 * (dm[i, k] + dm[j, k] - dm[i, j])

    if disallow_negative_branch_length and k_to_u < 0:
        k_to_u = 0

    return k_to_u


def _pair_members_to_new_node(dm, i, j, disallow_negative_branch_length):
    """Return the distance between a new node and decendants of that new node.
    Parameters
    ----------
    dm : skbio.DistanceMatrix
        The input distance matrix.
    i, j : str
        Identifiers of entries in the distance matrix to be collapsed (i.e.,
        the descendents of the new node, which is internally represented as
        `u`).
    disallow_negative_branch_length : bool
        Neighbor joining can result in negative branch lengths, which don't
        make sense in an evolutionary context. If `True`, negative branch
        lengths will be returned as zero, a common strategy for handling this
        issue that was proposed by the original developers of the algorithm.
    """
    n = dm.shape[0]
    i_to_j = dm[i, j]
    i_to_u = (0.5 * i_to_j) + ((dm[i].sum() - dm[j].sum()) / (2 * (n - 2)))

    if disallow_negative_branch_length and i_to_u < 0:
        i_to_u = 0

    j_to_u = i_to_j - i_to_u

    if disallow_negative_branch_length and j_to_u < 0:
        j_to_u = 0

    return i_to_u, j_to_u
    
    
