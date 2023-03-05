
mask = []
def square_sparse(A, ia, ja):
    for i in range(len(A)):
        for j in range(ia[i], ia[i+1]+1):
            friend = ja[j]
            k = friend
            while k < ia[friend+1]:
                friend_of_friend = ja[k]
                mask.append((i, friend_of_friend))
                k += 1
    