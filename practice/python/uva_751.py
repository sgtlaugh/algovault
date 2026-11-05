"""
Minimax with alpha beta pruning

AbstractNode: subclass with the game rules, get_utility() returns the value of a node under optimal play

Complexity: O(b^d) get_utility calls in the worst case, b = branching factor, d = depth of the game tree
Pruning only skips calls, how many depends on the order get_next_states returns, and without memoization
a state reachable through different move orders is searched again
Recursion depth d, each level keeps its list of up to b next states, O(b d) memory
"""

import math
from abc import ABC, abstractmethod


class AbstractNode(ABC):
    def __init__(self, state):
        self.state = state

    @abstractmethod
    def get_next_states(self):
        """
        :return: list of all next states reachable from the current state
        """
        pass

    @abstractmethod
    def is_terminal_node(self):
        """
        :return: True if the current node is terminal (i.e, game over) and False otherwise
        """
        pass

    @abstractmethod
    def get_terminal_node_utility(self):
        """
        :return: utility of the terminal node, for example, could be +1, -1, or 0 in the case of tic-tac-toe
        """
        pass

    def get_utility(self, is_first_player_turn=True, alpha=-math.inf, beta=math.inf):
        """
        Uses minimax with alpha beta pruning to get utility of any node

        :param is_first_player_turn: True if first player is to move, False otherwise
        :param alpha: alpha value (first player)
        :param beta: beta value (second player)
        :return: utility value of the node
        """
        if self.is_terminal_node():
            return self.get_terminal_node_utility()

        for next_state in self.get_next_states():
            value = next_state.get_utility(is_first_player_turn ^ True, alpha, beta)

            if is_first_player_turn and value > alpha:
                alpha = value

            if not is_first_player_turn and value < beta:
                beta = value

            if alpha >= beta:
                break

        return alpha if is_first_player_turn else beta


class NimNode(AbstractNode):
    """
    Pile of stones, a move takes 1 to 3 of them and whoever takes the last stone wins
    State is (stones left, whether the first player is to move)
    """

    def get_next_states(self):
        stones, first_to_move = self.state
        return [NimNode((stones - take, not first_to_move)) for take in range(1, min(3, stones) + 1)]

    def is_terminal_node(self):
        return self.state[0] == 0

    def get_terminal_node_utility(self):
        return -1 if self.state[1] else 1  # The player to move has no stones left, so the other player took the last one


import sys

EDGES = [(1, 2), (1, 3), (2, 3), (2, 4), (2, 5), (3, 5), (3, 6), (4, 5), (5, 6),
         (4, 7), (4, 8), (5, 8), (5, 9), (6, 9), (6, 10), (7, 8), (8, 9), (9, 10)]
BIT = {edge: 1 << i for i, edge in enumerate(EDGES)}
TRIANGLES = [BIT[(a, b)] | BIT[(a, c)] | BIT[(b, c)]
             for a, b, c in [(1, 2, 3), (2, 4, 5), (2, 3, 5), (3, 5, 6), (4, 7, 8), (4, 5, 8), (5, 8, 9), (5, 6, 9), (6, 9, 10)]]
FULL = (1 << len(EDGES)) - 1
SIDE_OF = {bit: [t for t in TRIANGLES if t & bit] for bit in BIT.values()}


def fill(state, bit):
    mask, a_turn, a, b = state
    mask |= bit
    gained = sum(1 for t in SIDE_OF[bit] if mask & t == t)
    if a_turn:
        a += gained
    else:
        b += gained
    return mask, a_turn if gained else not a_turn, a, b


class TriangleNode(AbstractNode):
    """
    State is (mask of filled lines, whether A is to move, triangles of A, triangles of B)
    """

    cache = {}

    def get_next_states(self):
        # A generator, so a cutoff skips building the remaining children
        mask = self.state[0]
        return (TriangleNode(fill(self.state, bit)) for bit in BIT.values() if not mask & bit)

    def is_terminal_node(self):
        return self.state[0] == FULL or max(self.state[2], self.state[3]) > len(TRIANGLES) // 2

    def get_terminal_node_utility(self):
        return 1 if self.state[2] > self.state[3] else -1

    def get_utility(self, is_first_player_turn=True, alpha=-1, beta=1):
        # A completed triangle gives another turn, so the mover comes from the state, not from alternating plies
        # Searched only with the window (-1, 1), which no cutoff ever narrows since values are -1 or 1,
        # so the result depends on the state alone and can be cached
        if self.state not in self.cache:
            self.cache[self.state] = super().get_utility(self.state[1], alpha, beta)
        return self.cache[self.state]


def main():
    tokens = sys.stdin.read().split()
    out, pos = [], 1
    for game in range(1, int(tokens[0]) + 1):
        m = int(tokens[pos])
        state = (0, True, 0, 0)
        for i in range(m):
            state = fill(state, BIT[(int(tokens[pos + 1 + 2 * i]), int(tokens[pos + 2 + 2 * i]))])
        pos += 1 + 2 * m

        res = TriangleNode(state).get_utility()
        out.append('Game %d: %s wins.' % (game, 'A' if res == 1 else 'B'))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
