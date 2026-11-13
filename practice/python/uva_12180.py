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


class GameNode(AbstractNode):
    """
    State is (mugs, diff, player, locks, swaps), seen from the player to move: mugs lists his mugs, then the
    opponent's, swaps his swaps left, then the opponent's, locks[j] the turns mug pair j stays unswappable
    diff is my score minus the opponent's so far, player 0 is me and maximizes
    A child is the board after a whole turn, extra moves included, so the players alternate as the template expects
    """

    m = 0
    table, turns = {}, {}

    def get_utility(self, is_first_player_turn=True, alpha=-math.inf, beta=math.inf):
        """
        Transposition table of bounds on the mover's future gain, shared by both players and every test case
        A position not in it is bounded by the stones left on the board, which halves the positions searched
        """
        mugs, diff, player, locks, swaps = self.state
        key = (mugs, locks, swaps)
        stones = sum(mugs)
        low, high = self.table.get(key, (-stones, stones))
        if player:
            low, high = -high, -low
        low, high = low + diff, high + diff

        if low == high or low >= beta:
            return low
        if high <= alpha:
            return high

        alpha, beta = max(alpha, low), min(beta, high)
        value = super().get_utility(is_first_player_turn, alpha, beta)
        if value <= alpha:
            high = value
        elif value >= beta:
            low = value
        else:
            low = high = value

        low, high = low - diff, high - diff
        self.table[key] = (-high, -low) if player else (low, high)
        return value

    def get_next_states(self):
        mugs, diff, player, locks, swaps = self.state
        key = (mugs, locks, swaps)
        if key not in self.turns:
            self.turns[key] = self.get_turns()

        sign = -1 if player else 1
        return [GameNode((child, diff + sign * gain, player ^ 1, new_locks, new_swaps))
                for gain, child, new_locks, new_swaps in self.turns[key]]

    def get_turns(self):
        """
        Every end of the mover's turn as (gain, mugs, locks, swaps) for the opponent, best gain first
        The mover's bowl sits at m and the opponent's at 2m + 1, both start at 0 so they hold the gains
        """
        mugs, diff, player, locks, swaps = self.state
        m, size = self.m, 2 * self.m + 2
        next_locks = tuple(x - 1 if x else 0 for x in locks)
        best = {}

        def finish(cells, new_locks, new_swaps):
            state = (cells[m + 1:-1] + cells[:m], new_locks, new_swaps)
            best[state] = max(best.get(state, -math.inf), cells[m] - cells[-1])

        if not any(mugs[:m]):
            finish(mugs[:m] + (0,) + mugs[m:] + (0,), next_locks, swaps[::-1])

        stack = [mugs[:m] + (0,) + mugs[m:] + (0,)]
        while stack:
            board = stack.pop()
            for i in range(m):
                stones = board[i]
                if not stones:
                    continue

                cur = list(board)
                cur[i] = 0
                last = (i + stones) % size
                was_empty = board[last] == 0
                for k in range(i + 1, i + stones + 1):
                    cur[k % size] += 1

                if last == m and any(cur[:m]):
                    stack.append(tuple(cur))
                    continue

                if last < m and was_empty and cur[2 * m - last]:
                    cur[m] += cur[last] + cur[2 * m - last]
                    cur[last] = cur[2 * m - last] = 0
                elif m < last < size - 1:
                    j = last - m - 1
                    if swaps[0] and not locks[j] and cur[j]:
                        swapped = list(cur)
                        swapped[j], swapped[last] = cur[last], cur[j]
                        finish(tuple(swapped), next_locks[:j] + (4,) + next_locks[j + 1:], (swaps[1], swaps[0] - 1))

                finish(tuple(cur), next_locks, swaps[::-1])

        return sorted(((gain, state[0], state[1], state[2]) for state, gain in best.items()), reverse=True)

    def is_terminal_node(self):
        return not any(self.state[0])

    def get_terminal_node_utility(self):
        return self.state[1]


def solve(m, cells):
    """
    Null windows closing in on the value reuse the table and visit 20% fewer positions than one wide window
    """
    GameNode.m = m
    mugs = tuple(cells[m + 1:-1] + cells[:m])
    diff, stones = cells[m] - cells[-1], sum(mugs)
    root = GameNode((mugs, diff, 1, (0,) * m, (3, 3)))

    low, high, guess = diff - stones, diff + stones, diff
    while low < high:
        beta = max(low + 1, min(high, guess))
        guess = root.get_utility(False, beta - 1, beta)
        if guess >= beta:
            low = guess
        else:
            high = guess
    return low


def main():
    tokens = sys.stdin.read().split()
    pos, out = 0, []
    while pos < len(tokens):
        t, pos = int(tokens[pos]), pos + 1
        for _ in range(t):
            m = int(tokens[pos])
            cells = [int(x) for x in tokens[pos + 1:pos + 3 + 2 * m]]
            pos += 3 + 2 * m
            out.append(str(solve(m, cells)))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
