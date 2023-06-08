"""
Utility functions

"""
# Copyright 2018 European Union
# This file is part of pyposeidon.
# Licensed under the EUPL, Version 1.2 or – as soon they will be approved by the European Commission - subsequent versions of the EUPL (the "Licence").
# Unless required by applicable law or agreed to in writing, software distributed under the Licence is distributed on an "AS IS" basis, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the Licence for the specific language governing permissions and limitations under the Licence.


from collections import deque


class Solution(object):
    def append_if(self, queue, x, y, island_counter):
        b = [[0, 1], [1, 0]], [[0, 1], [-1, 0]], [[0, -1], [1, 0]], [[0, -1], [-1, 0]]

        """Append to the queue only if in bounds of the grid and the cell value is 1."""
        if 0 <= x < len(self.grid) and 0 <= y < len(self.grid[0]):
            if self.grid[x][y] == "1":
                #                queue.append((x, y))
                for [i, j], [k, l] in b:
                    try:
                        if self.grid[x + i][y + j] in [
                            "1",
                            island_counter,
                        ] and self.grid[
                            x + k
                        ][y + l] in [
                            "1",
                            island_counter,
                        ]:
                            queue.append((x, y))
                            break
                    except:
                        pass

    def mark_neighbors(self, row, col, island_counter):
        """Mark all the cells in the current island with value = 2. Breadth-first search."""
        queue = deque()

        queue.append((row, col))
        while queue:
            x, y = queue.pop()
            self.grid[x][y] = island_counter

            self.append_if(queue, x - 1, y, island_counter)
            self.append_if(queue, x, y - 1, island_counter)
            self.append_if(queue, x + 1, y, island_counter)
            self.append_if(queue, x, y + 1, island_counter)
            self.append_if(queue, x - 1, y - 1, island_counter)
            self.append_if(queue, x + 1, y - 1, island_counter)
            self.append_if(queue, x + 1, y + 1, island_counter)
            self.append_if(queue, x - 1, y + 1, island_counter)

    def numIslands(self, grid):
        """
        :type grid: List[List[str]]
        :rtype: int
        """

        if not grid or len(grid) == 0 or len(grid[0]) == 0:
            return 0

        self.grid = grid

        row_length = len(grid)
        col_length = len(grid[0])

        island_counter = 0
        for row in range(row_length):
            for col in range(col_length):
                if self.grid[row][col] == "1":
                    # found an island
                    island_counter += 1

                    self.mark_neighbors(row, col, island_counter)

        return island_counter
