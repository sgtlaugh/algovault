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

LINES = [[4 * r + c for c in range(4)] for r in range(4)] + [[4 * r + c for r in range(4)] for c in range(4)]
LINES += [[0, 5, 10, 15], [3, 6, 9, 12]]


class TicTacToeNode(AbstractNode):
    """
    State is the 16 cells row by row, x always moves when both players have made the same number of moves
    """

    cache = {}

    def winner(self):
        for line in LINES:
            first = self.state[line[0]]
            if first != '.' and all(self.state[i] == first for i in line):
                return first
        return None

    def get_next_states(self):
        piece = 'x' if self.state.count('x') == self.state.count('o') else 'o'
        return [TicTacToeNode(self.state[:i] + piece + self.state[i + 1:]) for i in range(16) if self.state[i] == '.']

    def is_terminal_node(self):
        return self.winner() is not None or '.' not in self.state

    def get_terminal_node_utility(self):
        return {'x': 1, 'o': -1, None: 0}[self.winner()]

    def get_utility(self, is_first_player_turn=True, alpha=-math.inf, beta=math.inf):
        # Searched only with the window (0, 1), which no cutoff ever narrows since values are -1, 0 or 1,
        # so the result depends on the board alone and can be cached
        if self.state not in self.cache:
            self.cache[self.state] = super().get_utility(is_first_player_turn, alpha, beta)
        return self.cache[self.state]


def main():
    tokens = sys.stdin.read().split()
    out, pos = [], 0
    while tokens[pos] != '$':
        board = ''.join(tokens[pos + 1:pos + 5])
        pos += 5

        res = '#####'
        for i in range(16):
            if board[i] != '.':
                continue
            child = TicTacToeNode(board[:i] + 'x' + board[i + 1:])
            if child.get_utility(is_first_player_turn=False, alpha=0, beta=1) == 1:
                res = '(%d,%d)' % (i // 4, i % 4)
                break
        out.append(res)
    print('\n'.join(out))


if __name__ == '__main__':
    main()
