from collections import defaultdict
from itertools import combinations

import numpy as np
from scipy.special import comb


def compute_group_interaction_times(H_interaction, H_collaboration, size):
    interaction_groups = H_interaction.edges.members(dtype=dict)

    collaboration_groups = {
        frozenset(e) for e in H_collaboration.edges.filterby("size", size).members()
    }
    collaboration_non_groups = {
        frozenset(g) for g in combinations(list(H_interaction.nodes), size)
    }.difference(collaboration_groups)

    interaction_times = H_interaction.edges.attrs("duration", missing=0).asdict()

    t_c = 0
    t_nc = 0

    for cg in collaboration_groups:
        for idx, ig in interaction_groups.items():
            if cg.issubset(ig):
                t_c += interaction_times[idx] / (
                    len(ig) * len(collaboration_groups)
                )  # add group size and time

    for ncg in collaboration_non_groups:
        for idx, ig in interaction_groups.items():
            if ncg.issubset(ig):
                t_nc += interaction_times[idx] / (
                    len(ig) * len(collaboration_non_groups)
                )

    return t_c, t_nc


def compute_subgroup_interaction_times(
    H_interaction, H_collaboration, size, return_aggregates=False
):

    nodes = set(H_interaction.nodes).union(H_collaboration.nodes)

    interaction_groups = H_interaction.edges.members(dtype=dict)

    collaboration_groups = {
        frozenset(e) for e in H_collaboration.edges.filterby("size", size).members()
    }
    collaboration_non_groups = {
        frozenset(g) for g in combinations(list(nodes), size)
    }.difference(collaboration_groups)

    durations = H_interaction.edges.attrs("duration", missing=0).asdict()

    t_c = _count_subinteractions(
        interaction_groups,
        collaboration_groups,
        durations,
    )
    t_nc = _count_subinteractions(
        interaction_groups,
        collaboration_non_groups,
        durations,
    )

    if return_aggregates:
        a_c = _interpret_decomp1(t_c, size)
        a_nc = _interpret_decomp1(t_nc, size)
        return t_c, t_nc, a_c, a_nc
    else:
        return t_c, t_nc


def _count_subinteractions(groups1, groups2, durations):
    """This function adds up the time that subsets of a
    collaboration group spent in interaction groups.

    Parameters
    ----------
    groups1 : set of frozensets
        The interaction groups
    groups2 : set of frozensets
        The collaboration groups
    durations : dict of floats
        A dictionary where the keys are frozensets of the groups
        and the values are the time (normalized by the group size)

    Returns
    -------
    dict of dicts of floats
        The keys are the collaboration groups and the values are
        dictionaries where keys are the subgroups and values are the
        time spend interacting
    """
    interaction_dict = {}

    # This function looks for the biggest subset of the collaboration group
    # in an interaction group.
    def _iterate_over_subgroups(g1, g2):
        for size in range(len(g2), 1, -1):
            for sg in combinations(g2, size):
                if set(sg).issubset(g1):
                    return frozenset(sg)
        raise Exception("No subgroups!")

    for g2 in groups2:

        # the if statement is so that the times are
        # not reset to zero if there is a duplicated collaboration
        # group (should never happen though!)
        if g2 not in interaction_dict:
            interaction_dict[g2] = {}
            for size in range(2, len(g2) + 1):
                for sg in combinations(g2, size):
                    interaction_dict[g2][frozenset(sg)] = 0

        # count up subsets of the collab groups that are in the
        # interaction groups.
        for idx, g1 in groups1.items():
            try:
                sg = _iterate_over_subgroups(g1, g2)
                interaction_dict[g2][sg] += durations[idx] / len(g1)
            except:
                pass
    return interaction_dict


def _interpret_decomp1(d, size):
    interaction_types = defaultdict(lambda: 0)
    # iterate over all collaboration groups and
    # look at the subinteractions
    n = len([g for g in d if len(g) == size])
    print(n)
    for g, val in d.items():
        # only look at the particular sized collaboration
        # group
        if len(g) == size:
            # iterate over subgroup sizes from 2 to the size of the collab group
            for s in range(2, size + 1):
                # get the list of interacting subgroups of a certain size
                a = np.array(
                    [
                        val[frozenset(sg)]
                        for sg in combinations(g, s)
                        if val[frozenset(sg)] != 0
                    ]
                )
                # count them
                n_sgi = np.count_nonzero(a)
                # this loop iteratively finds the smallest common time,
                # e.g., if (1, 2), (1, 3), and (2, 3) interacted for 0.25, 0.5, and 0.75 hrs,
                # then the smallest common time is 0.25
                # and assigns this time to a collection of iteractions where 3 groups of two interact
                # we then subtract off that time, remove the zero and do this iteratively, until no time
                # is left.
                while n_sgi > 0 and len(a):
                    t = a.min()
                    a -= t
                    interaction_types[(s, n_sgi)] += float(t) / n
                    a = a[np.nonzero(a)]
                    sgs = np.count_nonzero(a)

    return dict(interaction_types)


def create_data_matrix(t_nc, t_c, size):
    interaction_types = defaultdict(lambda: 0)

    n_c = len([g for g in t_c if len(g) == size])
    n_nc = len([g for g in t_nc if len(g) == size])

    idx = {}
    it = 0
    for s in range(size, 1, -1):
        for i in range(comb(size, s, exact=True)):
            idx[(s, i + 1)] = it
            it += 1

    A = np.zeros((n_c + n_nc, len(idx) + 1))

    # iterate over all collaboration groups and
    # look at the subinteractions
    for i, (g, val) in enumerate(t_nc.items()):
        # only look at the particular sized collaboration
        # group
        if len(g) == size:
            # iterate over subgroup sizes from 2 to the size of the collab group
            for s in range(2, size + 1):
                # get the list of interacting subgroups of a certain size
                a = np.array(
                    [
                        val[frozenset(sg)]
                        for sg in combinations(g, s)
                        if val[frozenset(sg)] != 0
                    ]
                )
                # count them
                nz = np.count_nonzero(a)
                # this loop iteratively finds the smallest common time,
                # e.g., if (1, 2), (1, 3), and (2, 3) interacted for 0.25, 0.5, and 0.75 hrs,
                # then the smallest common time is 0.25
                # and assigns this time to a collection of iteractions where 3 groups of two interact
                # we then subtract off that time, remove the zero and do this iteratively, until no time
                # is left.
                while nz > 0:
                    t = a.min()
                    a -= t
                    j = idx[(s, nz)]
                    A[i, j] = float(t)
                    a = a[np.nonzero(a)]
                    nz = np.count_nonzero(a)
                A[i, -1] = 0

    for i, (g, val) in enumerate(t_c.items()):
        # only look at the particular sized collaboration
        # group
        if len(g) == size:
            # iterate over subgroup sizes from 2 to the size of the collab group
            for s in range(2, size + 1):
                # get the list of interacting subgroups of a certain size
                a = np.array(
                    [
                        val[frozenset(sg)]
                        for sg in combinations(g, s)
                        if val[frozenset(sg)] != 0
                    ]
                )
                # count them
                nz = np.count_nonzero(a)

                # this loop iteratively finds the smallest common time,
                # e.g., if (1, 2), (1, 3), and (2, 3) interacted for 0.25, 0.5, and 0.75 hrs,
                # then the smallest common time is 0.25
                # and assigns this time to a collection of iteractions where 3 groups of two interact
                # we then subtract off that time, remove the zero and do this iteratively, until no time
                # is left.
                while nz > 0:
                    t = a.min()
                    a -= t
                    j = idx[(s, nz)]
                    A[i, j] = float(t)
                    a = a[np.nonzero(a)]
                    nz = np.count_nonzero(a)
                A[i, -1] = 1

    return A, idx
