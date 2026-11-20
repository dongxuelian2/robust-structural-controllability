"""Vertex Cover reductions with a backup shared across a matching."""
from dataclasses import dataclass
from itertools import combinations

from assignment_transfer import branching_presentation, normalize_dominated
from quota_cotransversal import quota_cotransversal_graph, quota_from_placement


@dataclass(frozen=True)
class BackupCover:
    vertices: int
    jobs: int
    options: tuple
    forced: frozenset
    groups: tuple


def shared_backup_instance(vertices, edges, packed, triangles=(), strong=False):
    """Build an assignment with one job per edge outside the triangles.

    Edges use distinct endpoints in range(vertices). Packed edges form a
    matching disjoint from the triangle edges. Triangles are vertex-disjoint.
    Each vertex must occur in a retained edge or in a strong triangle job.
    """
    edges = {tuple(sorted(edge)) for edge in edges}
    packed = tuple(tuple(sorted(edge)) for edge in packed)
    triangles = tuple(tuple(triangle) for triangle in triangles)
    triangle_edges = {tuple(sorted(edge)) for triangle in triangles
                      for edge in combinations(triangle, 2)}
    packed_ends = [v for edge in packed for v in edge]
    triangle_ends = [v for triangle in triangles for v in triangle]
    if (any(len(edge) != 2 or not 0 <= edge[0] < edge[1] < vertices for edge in edges)
            or not packed or not set(packed) <= edges - triangle_edges
            or len(packed_ends) != len(set(packed_ends))
            or any(len(triangle) != 3 for triangle in triangles)
            or len(triangle_ends) != len(set(triangle_ends))
            or not triangle_edges <= edges):
        raise ValueError('Use a matching and vertex-disjoint triangles from the source graph.')
    remaining = tuple(sorted(edges - triangle_edges - set(packed)))
    return guarded_backup_assignment(vertices, remaining, packed, triangles, strong)


def guarded_backup_assignment(vertices, private_jobs, packed_jobs, triangles=(), strong=False,
                              *, optional=()):
    """Build the shared backup and private guards. Optional indices need no forcing."""
    optional = frozenset(optional)
    if not optional <= set(range(len(private_jobs))):
        raise ValueError('Optional guards must index private jobs.')
    neighborhoods = private_jobs + packed_jobs + (triangles if strong else ())
    guard = len(neighborhoods)
    options = [0] * vertices
    for j, neighbors in enumerate(neighborhoods):
        for v in neighbors:
            options[v] |= 1 << j
    if not all(options):
        raise ValueError('Every vertex must serve at least one retained job.')

    # The guard forces both the shared backup and its guard-only partner.
    backup = sum(1 << j for j in range(len(private_jobs), len(private_jobs) + len(packed_jobs)))
    options.extend((backup | (1 << guard), 1 << guard))
    forced = {vertices, vertices + 1}
    jobs = guard + 1
    if strong:
        for j in range(len(private_jobs)):
            if j in optional:
                options.append(1 << j)
                continue
            forced.update((len(options), len(options) + 1))
            private_guard = 1 << jobs
            options.extend((private_guard | (1 << j), private_guard))
            jobs += 1
        groups = ()
    else:
        options.extend(1 << j for j in range(len(private_jobs)))
        private = [vertices + 2 + j for j in range(len(private_jobs)) if j not in optional]
        forced.update(private)
        if len(private) % 2:
            private.append(vertices + 1)
        else:
            private.extend((vertices, vertices + 1))
        groups = triangles + tuple(zip(private[::2], private[1::2]))
    return BackupCover(vertices, jobs, tuple(options), frozenset(forced), groups)


def shared_backup_quota_graph(source, length=1):
    presentation, relaxed = branching_presentation(source.jobs, source.options)
    return quota_cotransversal_graph(presentation, len(relaxed), source.groups, length)


def cover_from_assignment(source, selected):
    """Extract the cover after original-resource normalization."""
    assert source.forced <= selected
    return frozenset(v for v in selected if v < source.vertices)


def cover_from_quota_placement(source, target, selected):
    support = quota_from_placement(target, selected)
    _, relaxed = branching_presentation(source.jobs, source.options)
    original = normalize_dominated(source.jobs, source.options, relaxed, support)
    assert all(len(original.intersection(group)) >= 2 for group in source.groups)
    return cover_from_assignment(source, original)
