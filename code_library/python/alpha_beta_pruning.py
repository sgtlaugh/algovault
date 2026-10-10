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
        super().__init__()

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

            if is_first_player_turn and value > alpha:  # First player maximizes score
                alpha = value

            if not is_first_player_turn and value < beta:  # Second player minimizes score
                beta = value

            if alpha >= beta:  # Alpha beta pruning applied here
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


def main():
    # The first player loses Nim with moves of 1 to 3 stones exactly when the pile is a multiple of 4
    for stones in range(1, 16):
        assert NimNode((stones, True)).get_utility(is_first_player_turn=True) == (-1 if stones % 4 == 0 else 1)


if __name__ == '__main__':
    main()
