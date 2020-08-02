import math
from abc import ABC, abstractmethod


class AbstractNode(ABC):
    @abstractmethod
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


class Node(AbstractNode):
    def __init__(self, state):
        super().__init__(state)

    def get_next_states(self):
        pass

    def is_terminal_node(self):
        pass

    def get_terminal_node_utility(self):
        pass


def main():
    state = {}
    start_node = Node(state)

    # Throws an exception because the Node methods are not implemented
    node_value = start_node.get_utility(is_first_player_turn=True)


if __name__ == '__main__':
    main()
