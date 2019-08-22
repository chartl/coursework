"""\

Pileup utilities

"""
from collections import Counter

_CIG_DELETION, _CIG_INSERTION = 2, 1


def next_read(read_seq):
    try:
        return read_seq.next()
    except StopIteration:
        return None


def make_pileups(read_seq, reference, contig, offset, window_size=0, stop=None):
    """\
    Generate a sequence of Pileup objects from a given sequence of
    ***sorted*** reads, starting at the given reference offset.

    """
    if offset - window_size < 0:
        raise ValueError('Cannot use a reference window that goes off-reference')
    read_seq = iter(read_seq)  # just in case
    # get to the first read
    cur_read = next_read(read_seq)
    while cur_read and cur_read.aend < offset:  # 0 + 50 = 50; offset is 49
        cur_read = read_seq.next()  # this is why bams have indexes...
    # now get all the reads that overlap the position, and put into a list
    alive_elements = []
    stop = stop or len(reference)
    while offset < stop:
        alive_elements = [e for e in alive_elements if not e.is_last_position()]
        for e in alive_elements:
            e._next_position()
        new_reads = []
        while cur_read and cur_read.pos <= offset:
            if not cur_read.is_unmapped:
                new_reads.append(cur_read)  # assumes reads are the same length
            cur_read = next_read(read_seq)
        alive_elements.extend([PileupElement(read, offset-read.pos) for read in new_reads])
        # make the reference window
        win_size = min(window_size, len(reference) - offset)
        ref_win = reference[offset - win_size:offset + win_size + 1]
        yield Pileup(contig, offset, ref_win, alive_elements)
        offset += 1


class Pileup(object):
    """\
    Information about the read sequences at a particular locus. This is a wrapper
    around a sequence of pileup elements, plus some extra fluff about reference
    location, etc.

    """
    def __init__(self, contig, position, ref_window, elements):
        self.contig = contig 
        self.position = position
        self.reference_window = ref_window
        self.reference_base = self.reference_window[len(self.reference_window)/2]
        # AACTT
        #   ^
        self.elements = list(elements)

    def __str__(self):
        return '(Pileup pos={} ref={}) {}'.format(self.position, self.reference_base,
                                                  ''.join([e.base() for e in self.elements if e.base() is not None]))

    def __iter__(self):
        return iter(self.elements)


    def __len__(self):
        return len(self.elements)

    @property
    def base_counts(self):
        return Counter((e.base() for e in self.elements))

    @property
    def counts(self):
        return self.base_counts 


class PileupElement(object):
    """\
    Class encapsulating ``vertical'' information at a locus, but providing access
    to the underlying read

    TODO: get insertion bases
    TODO: get deletion size

    """
    def __init__(self, read, offset):
        if read.cigar[0][0] in {_CIG_DELETION, _CIG_INSERTION}:
            raise ValueError('Reads cannot start with insertion or deletion')
        if read.cigar[-1][0] in {_CIG_DELETION, _CIG_INSERTION}:
            raise ValueError('Reads cannot end with insertion or deletion (seriously?)')
        self.underlying_read = read
        if offset == 0:
            self.read_offset = offset
            self._cig_offset, self._cig_elem_offset = 0, 0
        else:
            try:
                self.read_offset, self._cig_offset, \
                  self._cig_elem_offset = self._advance_to_offset(read, offset)
            except ValueError:
                raise ValueError('Error attempting to advance read\n{}\n'
                                 'at pos {} to offset {}'.format(str(read), read.pos, offset))

    def _advance_to_offset(self, read, ref_base_offset):
        """\
        Creating an element from the middle of a read. This means keeping track of
        where the sequence elements line up with regards to the cigar and all that.

        For now, and for safety, start the read off at 0 and advance `manually' to
        the desired offset.

        """
        self.read_offset = 0
        self._cig_elem_offset, self._cig_offset = 0, 0
        for _ in xrange(ref_base_offset):
            self._next_position()
        return self.read_offset, self._cig_offset, self._cig_elem_offset

    def _cigarelement(self):
        try:
            return self.underlying_read.cigar[self._cig_offset]
        except IndexError:
            raise ValueError('Error in pileup: attempted to get index {} of cigar {} for read {}'.format(self._cig_offset, self.underlying_read.cigar, self.underlying_read)) 

    def base(self):
        if self._cigarelement()[0] in {_CIG_DELETION, _CIG_INSERTION}:
            return None
        return self.underlying_read.seq[self.read_offset]

    def base_qual(self):
        if self._cigarelement()[0] in {_CIG_DELETION, _CIG_INSERTION}:
            return None
        return self.underlying_read.qual[self.read_offset]

    def deletion_next(self):
        if self._cig_elem_offset == len(self.underlying_read.cigar)-1:
            # if we're the last element, next one can't be a deletion
            return False
        if self._cig_elem_offset == self.underlying_read.cigar[self._cig_offset][1]-1:
            # if at the very last position of a cigar element
            if self.underlying_read.cigar[self._cig_offset+1][0] == _CIG_DELETION:
                # if the next element is a deletion
                    return True
        return False

    def insertion_next(self):
        if self._cig_elem_offset == len(self.underlying_read.cigar)-1:
            # if we're the last element, next one can't be an insertion
            return False
        if self._cig_elem_offset == self.underlying_read.cigar[self._cig_offset][1]-1:
            # if at the very last position of a cigar element
            if self.underlying_read.cigar[self._cig_offset+1][0] == _CIG_INSERTION:
                # if the next element is a deletion
                    return True
        return False

    def is_last_position(self):
        """\
        Whether or not this is the last aligned base in the read.

        Note that this returns TRUE if the only remaining cigar element is an
        insertion

        """
        if self.read_offset == len(self.underlying_read.seq) - 1:
            return True
        if self.underlying_read.cigar[-1][1] == _CIG_INSERTION:
            # read ends with insertion (RARE)
            if self._cig_offset == len(self.underlying_read.cigar) - 2:
                # at the 2nd-to-last cigar element
                return self._cig_elem_offset == self._cigarelement()[1] - 1  # last base
        return False

    def _next_position(self):
        """\
        Advance to the next *reference* position.

        This means skipping any insertions!

        """
        if self.is_last_position():
            raise ValueError('Pileup element attempted to exceed read length:'
                             '\n{}'.format(str(self.underlying_read)))
        # first deal with the cigar elements
        if self._cig_elem_offset == self._cigarelement()[1] - 1:
            # we are at the end of the cigar, so need to increment
            if self._cigarelement()[0] != _CIG_DELETION:
                self.read_offset += 1
            self._cig_offset += 1
            self._cig_elem_offset = 0
            if self._cigarelement()[0] == _CIG_INSERTION:
                self.read_offset += self._cigarelement()[1]
                self._cig_offset += 1
        # update the read offset if necessary
        else:
            self._cig_elem_offset += 1
            if self._cigarelement()[0] != _CIG_DELETION:
                self.read_offset += 1
