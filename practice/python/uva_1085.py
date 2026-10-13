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

# Slots 0..7 hold the peaks, a floor on DOWN[k][0] rests on peaks DOWN[k][1], DOWN[k][2] and a peak pair UP[k][0], UP[k][1]
# rests on floor UP[k][2]
DOWN = [(10, 1, 2), (13, 3, 4), (16, 5, 6), (19, 9, 11), (22, 12, 14), (25, 18, 20)]
UP = [(i, i + 1, i + 2) for i in range(8, 24, 3)]


def triangle_score(a, b, c):
    res = abs(a) + abs(b) + abs(c)
    return -res if (a < 0) + (b < 0) + (c < 0) >= 2 else res


class HouseNode(AbstractNode):
    """
    State is (cards, held, turn, player, score): cards[slot] is a signed card (red > 0, black < 0, 0 for empty),
    held[p] the card player p holds, player 0 is Axel, score is Axel's minus Birgit's from the triangles built so far
    """

    deck = []

    def get_next_states(self):
        cards, held, turn, p, score = self.state
        card, mine, res = self.deck[turn], held[p], []
        nxt, q = turn + 1, p ^ 1

        def keep(value):
            return (value, held[1]) if p == 0 else (held[0], value)

        if not mine:
            res.append(HouseNode((cards, keep(card), nxt, q, score)))

        # A floor takes the drawn card, or the held one while the drawn card is kept
        for placed, kept in ((card, mine), (mine, card)) if mine else ((card, mine),):
            new_held = keep(kept)
            for floor, left, right in DOWN:
                if not cards[floor] and cards[left] and cards[right]:
                    gain = triangle_score(placed, cards[left], cards[right])
                    res.append(HouseNode((cards[:floor] + (placed,) + cards[floor + 1:], new_held, nxt, q, score + gain)))

        if mine:
            new_held = keep(0)
            for left, right, floor in UP:
                if cards[floor] and not cards[left] and not cards[right]:
                    gain = score + triangle_score(card, mine, cards[floor])
                    res.append(HouseNode((cards[:left] + (card, mine) + cards[right + 1:], new_held, nxt, q, gain)))
                    res.append(HouseNode((cards[:left] + (mine, card) + cards[right + 1:], new_held, nxt, q, gain)))

        # Trying the greediest move first tightens alpha and beta early, which prunes far more of the tree
        res.sort(key=lambda node: node.state[4], reverse=p == 0)
        return res

    def is_terminal_node(self):
        return self.state[2] == len(self.deck)

    def get_terminal_node_utility(self):
        held = self.state[1]
        return self.state[4] + held[0] + held[1]


def main():
    tokens = sys.stdin.read().split()
    pos, case, out = 0, 0, []
    while tokens[pos] != "End":
        name, m = tokens[pos], int(tokens[pos + 1])
        deck = [int(s[:-1]) * (1 if s[-1] == 'R' else -1) for s in tokens[pos + 2:pos + 2 + 2 * m]]
        pos += 2 + 2 * m

        HouseNode.deck = deck
        root = HouseNode((tuple(deck[:8]) + (0,) * 18, (0, 0), 8, 0 if deck[0] > 0 else 1, 0))
        diff = root.get_utility(is_first_player_turn=deck[0] > 0)

        case += 1
        if diff == 0:
            out.append("Case %d: Axel and Birgit tie" % case)
        else:
            gain = diff if name == "Axel" else -diff
            out.append("Case %d: %s %s %d" % (case, name, "wins" if gain > 0 else "loses", abs(gain)))
    print("\n".join(out))


if __name__ == '__main__':
    main()
