import math

from stress import rng, scaled

from alpha_beta_pruning import AbstractNode, NimNode


class TreeNode(AbstractNode):
    """State is a leaf utility or a list of child states, visits counts the nodes the search expands"""
    visits = 0

    def get_next_states(self):
        TreeNode.visits += 1
        return [TreeNode(child) for child in self.state]

    def is_terminal_node(self):
        return not isinstance(self.state, list)

    def get_terminal_node_utility(self):
        return self.state


def random_tree(depth):
    if depth == 0 or rng.randint(0, 5) == 0:
        return rng.randint(-20, 20)
    return [random_tree(depth - 1) for _ in range(rng.randint(1, 4))]


def minimax(state, maximizing):
    if not isinstance(state, list):
        return state
    values = [minimax(child, not maximizing) for child in state]
    return max(values) if maximizing else min(values)


def count_inner(state):
    return isinstance(state, list) and 1 + sum(count_inner(child) for child in state)


def main():
    pruned = 0
    for _ in range(scaled(600)):
        tree, first = random_tree(rng.randint(0, 7)), rng.randint(0, 1) == 1
        exact = minimax(tree, first)

        TreeNode.visits = 0
        assert TreeNode(tree).get_utility(first) == exact
        assert TreeNode.visits <= count_inner(tree)
        pruned += count_inner(tree) - TreeNode.visits

        # Fail-hard window semantics: exact inside (alpha, beta), otherwise on the correct side of the window
        alpha = rng.choice([-math.inf, rng.randint(-25, 25)])
        beta = rng.choice([math.inf, rng.randint(-25, 25)])
        if alpha < beta:
            value = TreeNode(tree).get_utility(first, alpha, beta)
            if alpha < exact < beta:
                assert value == exact
            elif exact <= alpha:
                assert value <= alpha
            else:
                assert value >= beta
    assert pruned > 0, 'alpha beta must skip some subtrees'

    for stones in range(16):
        for first in (True, False):
            loser = stones % 4 == 0  # the player to move loses on multiples of 4
            winner_is_first = first != loser
            assert NimNode((stones, first)).get_utility(first) == (1 if winner_is_first else -1)


if __name__ == '__main__':
    main()
