from collections import deque


class Solution(object):
    def append_if(self, queue, x, y):
        """Append to the queue only if in bounds of the grid and the cell value is 1."""
        if 0 <= x < len(self.grid) and 0 <= y < len(self.grid[0]):
            if self.grid[x][y] == '1':
                queue.append((x, y))

    def append_if_(self, queue, x, y, i):
        """Append to the queue only if in bounds of the grid and the cell value is 1."""
        if 0 <= x < len(self.grid) and 0 <= y < len(self.grid[0]):
            if self.grid[x][y] == '1':
              if 0 <= x+i[0] < len(self.grid) and 0 <= y+i[1] < len(self.grid[0]):
                 if self.grid[x+i[0]][y+i[1]] == '1' : 
                     queue.append((x, y))
              elif 0 <= x+i[2] < len(self.grid) and 0 <= y+i[3] < len(self.grid[0]):
                 if self.grid[x+i[2]][y+i[3]] == '1' : 
                     queue.append((x, y))
#             elif 0 <= x+i[4] < len(self.grid) and 0 <= y+i[5] < len(self.grid[0]):
#                if self.grid[x+i[4]][y+i[5]] == '1' :
#                    queue.append((x, y))

    def mark_neighbors(self, row, col, island_counter):
        """Mark all the cells in the current island with value = 2. Breadth-first search."""
        queue = deque()

        queue.append((row, col))
        while queue:
            x, y = queue.pop()
            self.grid[x][y] = island_counter

#           self.append_if_(queue, x - 1, y, [-1,0,0,1,0,-1])
#           self.append_if_(queue, x, y - 1, [0,-1,-1,0,1,0])
#           self.append_if_(queue, x + 1, y, [1,0,0,-1,0,1])
#           self.append_if_(queue, x, y + 1, [0,1,-1,0,1,0])

            self.append_if_(queue, x - 1, y, [0,1,0,-1])
            self.append_if_(queue, x, y - 1, [-1,0,1,0])
            self.append_if_(queue, x + 1, y, [0,1,0,-1])
            self.append_if_(queue, x, y + 1, [-1,0,1,0])
#            self.append_if(queue, x - 1, y - 1)
#            self.append_if(queue, x + 1, y - 1)
#            self.append_if(queue, x + 1, y + 1)
#            self.append_if(queue, x - 1, y + 1)


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
                if self.grid[row][col] == '1':
                    # found an island
                    island_counter += 1

                    self.mark_neighbors(row, col, island_counter)

        return island_counter
