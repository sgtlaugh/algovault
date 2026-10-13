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


import sys


class PawnNode(AbstractNode):
    """
    State is (board, white to move, plies played), board a string of 16 cells row by row from the top
    White pawns 'P' move up, black pawns 'p' move down
    """

    def get_next_states(self):
        """
        Promotions first, then captures, then pushes, most advanced pawns first within each group
        Without memoization the search lives on pruning: this order visits about 8x fewer nodes than board order
        """
        if hasattr(self, 'children'):
            return self.children

        board, white, depth = self.state
        me, enemy, step = ('P', 'p', -4) if white else ('p', 'P', 4)
        moves = []
        for i in (range(16) if white else range(15, -1, -1)):
            if board[i] != me or not 0 <= i + step < 16:
                continue
            for dc, target in ((0, '.'), (-1, enemy), (1, enemy)):
                j = i + step + dc
                if 0 <= i % 4 + dc < 4 and board[j] == target:
                    cells = list(board)
                    cells[i], cells[j] = '.', me
                    priority = 0 if not 0 <= j + step < 16 else (1 if dc else 2)
                    moves.append((priority, PawnNode((''.join(cells), not white, depth + 1))))

        moves.sort(key=lambda move: move[0])
        self.children = [node for _, node in moves]
        return self.children

    def is_terminal_node(self):
        board = self.state[0]
        return 'P' in board[:4] or 'p' in board[12:] or not self.get_next_states()

    def get_terminal_node_utility(self):
        """
        The game ends on the loser's turn: they were promoted against or are stalemated
        A faster win scores higher for white and a slower loss scores higher for black
        """
        board, white, depth = self.state
        return depth - 100 if white else 100 - depth


def main():
    tokens = sys.stdin.read().split()
    for case in range(int(tokens[0])):
        board = ''.join(tokens[1 + 4 * case:5 + 4 * case])
        value = PawnNode((board, True, 0)).get_utility(is_first_player_turn=True)
        moves = 100 - abs(value)
        print('%s (%d)' % ('white' if value > 0 else 'black', moves))


if __name__ == '__main__':
    main()
