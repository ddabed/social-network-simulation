import numpy as np 
import py2neo as pn 
import random
import math
import copy


def create_ba_subgraph(M, m_0, n, relation):
    # Generate Users-Friends un-directed subgraph
    # Generate Users-Follower directed subgraph
    node_list = []
    for i in range(m_0):
        node_list.append(pn.Node('User', name='user' + str(i)))
    rel_list = []
    for i in range(len(node_list)):
        for j in range((i+1),len(node_list)):
            rel_list.append(pn.Relationship(node_list[i], relation, node_list[j]))

    # 1. Creat node list by node function
    # 2. Creat relation list by Relation function
    # 3. Creat subgraph by Subgraph function
    t = M - m_0     # number of iteration
    k = []          # save nodes degree
    k_0 = m_0 - 1
    p_k = []        # save nodes priority probability
    p_0 = 1/m_0
    k_all = 0

    for i in range(m_0):
        p_k.append(p_0)
        k.append(k_0)
        k_all += k_0

    for i in range(t):
        m_0_t = m_0 + i     # number of nodes at time t
        m_0_1 = m_0 + i - 1 # number of nodes at time t-1
        node_list.append(pn.Node('User', name='user' + str(m_0_t)))  # add one node
        add_edge = 1
        j_choose = -1
        while(add_edge <= n):
            for j in range(m_0_t):
                if j != j_choose:   # to ensure no repeated edges
                    p_random = random.random()
                    if p_random <= p_k[j] and add_edge <= n:
                        j_choose = j
                        k_j = k[j]
                        p_k_j = p_k[j]
                        r_random = random.random()
                        if r_random > 0.5:
                            rel_list.append(pn.Relationship(node_list[j], relation, node_list[-1]))
                        else:
                            rel_list.append(pn.Relationship(node_list[-1], relation, node_list[j]))
                        add_edge += 1
                        k[j] = k_j + 1
                        k_all += 2
                        p_k[j] = (k_j + 1)/k_all
        k.append(n)
        p_k.append(n/k_all)
    s = pn.Subgraph(node_list,rel_list)
    return s

def create_pp_ba_subgraph(M, m_0, n, post_u0, N, graph):
    post_u0_list = np.random.choice(N,m_0)
    post_0_list = list(range(m_0))
    user_list = [pn.Node('User', name='user' + str(i)) for i in N]

    k = [0] * len(N)    # list save the number of posts by users
    node_list = []
    for i in post_0_list:
        node_list.append(pn.Node('Post', name='post' + str(i)))
    rel_list = []

    for i in range(len(node_list)):
        p_u = np.random.choice(post_u0_list,1)[0]
        k[p_u] += 1
        rel_list.append(pn.Relationship(user_list[p_u], 'published', node_list[i]))


    # 1. Creat node list by Node function
    # 2. Creat relation list by Relation function
    # 3. Creat subgraph by Subgraph function
    t = M - m_0     # number of iteration
    k_all_friends = 0
    for i in N:
        friends = list(graph.run('MATCH (n:User {name:"user' + str(i) + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
        k_all_friends += friends

    k_all_follow = 0
    for i in N:
        follow = list(graph.run('MATCH (n:User {name:"user' + str(i) + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
        k_all_follow += follow

    for i in range(t):
        m_0_t = m_0 + i     # number of nodes at time t
        node_list.append(pn.Node('Post', name='post' + str(m_0_t)))  # add one node
        p_j = [0] * len(N)  # save list of probability
        for j in N:
            p_j_friends = list(graph.run('MATCH (n:User {name:"user' + str(j) + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
            p_j_follow = list(graph.run('MATCH (n:User {name:"user' + str(j) + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
            p_j[j] = (p_j_friends + p_j_follow + k[j]) / (k_all_follow + k_all_friends + sum(k))
        user = np.random.choice(user_list,1, p=p_j)[0]   # roulette wheel selection
        rel_list.append(pn.Relationship(user, 'published', node_list[-1]))
        k[user_list.index(user)] += 1
    s = pn.Subgraph(node_list,rel_list)
    return s

def create_pv_ba_subgraph(P,N,graph):
    # depends on friends/followers/reads
    user_list = [pn.Node('User', name='user' + str(i)) for i in N]
    post_list = [pn.Node('Post', name='post' + str(i)) for i in P]
    view_list = [0] * len(P)
    rel_list = []
    k_all_friends = 0
    k_all_follow = 0
    for i in P:
        user = list(graph.run('match (n:Post{name:"post' + str(i) + '"})<-[:published]-(a) return a').data()[0].values())[0].nodes[0]['name']
        friends = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
        k_all_friends += friends
        follow = list(graph.run('MATCH (n:User {name:"' + user + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
        k_all_follow += follow
    for i in P:
        user = list(graph.run('match (n:Post{name:"post' + str(i) + '"})<-[:published]-(a) return a').data()[0].values())[0].nodes[0]['name']
        user_list_m = copy.deepcopy(user_list)
        for n in user_list_m:
            if n['name'] == user:
                user_list_m.remove(n)
        friends = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
        follow = list(graph.run('MATCH (n:User {name:"' + user + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
        p = (friends + follow + view_list[i]) / (k_all_friends + k_all_follow + sum(view_list))
        p_random = random.random()
        if p_random <= p:
            user_choice = np.random.choice(user_list_m,1)[0]
            rel_list.append(pn.Relationship(user_choice,'viewed',post_list[i]))
            view_list[i] += 1
    s = pn.Subgraph(post_list,rel_list)
    return s

def create_pl_ba_subgraph(P,N,graph):
    user_list = [pn.Node('User', name='user' + str(i)) for i in N]
    post_list = [pn.Node('Post', name='post' + str(i)) for i in P]
    liked_list = [0] * len(P)
    rel_list = []
    k_all_friends = 0
    k_all_follow = 0
    k_all_posted = 0
    for i in P:
        user = list(graph.run('match (n:Post{name:"post' + str(i) + '"})<-[:published]-(a) return a').data()[0].values())[0].nodes[0]['name']
        friends = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
        k_all_friends += friends
        follow = list(graph.run('MATCH (n:User {name:"' + user + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
        k_all_follow += follow
        post = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:published]->(a) return count(a)').data()[0].values())[0] 
        k_all_posted += post
    for i in P:
        user = list(graph.run('match (n:Post{name:"post' + str(i) + '"})<-[:published]-(a) return a').data()[0].values())[0].nodes[0]['name']
        user_list_m = copy.deepcopy(user_list)
        for n in user_list_m:
            if n['name'] == user:
                user_list_m.remove(n)
        friends = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:friends]-(a) return count(a)').data()[0].values())[0]
        follow = list(graph.run('MATCH (n:User {name:"' + user + '"})<-[:follow]-(a) return count(a)').data()[0].values())[0]
        post = list(graph.run('MATCH (n:User {name:"' + user + '"})-[:published]->(a) return count(a)').data()[0].values())[0]
        for j in N:
            p = (friends + follow + post + liked_list[i]) / (k_all_friends + k_all_follow + k_all_posted + sum(liked_list))
            p_random = random.random()
            if p_random <= p:
                user_choice = np.random.choice(user_list_m,1)[0]
                rel_list.append(pn.Relationship(user_choice,'liked',post_list[i]))
                liked_list[i] += 1
    s = pn.Subgraph(post_list,rel_list)
    return s



if __name__ == '__main__':
    users = 20      # users nodes
    posts = 30       # posts
    post_0, post_u0 = 3,3
    m_0 = 3         # initial nodes
    n = 2           # every time the new node connect to n known nodes, n<=n_user0
    N = list(range(users))
    P = list(range(posts))
    # start a new project in Neo4j and set connections
    graph = pn.Graph(
        host = 'localhost',
        http_port = '7474',
        user = 'neo4j',
        password = '2500'
        )    
    # stage 1
    s_user_friend = create_ba_subgraph(users, m_0, n, relation='friends')
    # stage 2
    s_user_follow = create_ba_subgraph(users, m_0, n, relation='follow')
    graph.run('match (n:User) detach delete n')
    graph.create(s_user_friend)
    # stage 3
    graph.merge(s_user_follow,'User','name')

    # stage 4
    s_posts_publish = create_pp_ba_subgraph(posts, post_0, n, post_u0, N, graph)
    # stage 5
    graph.merge(s_posts_publish,'User','name')
    # stage 6
    s_posts_viewed = create_pv_ba_subgraph(P,N,graph)
    # stage 7
    graph.merge(s_posts_viewed,'User','name')
    # stage 8
    s_posts_liked = create_pl_ba_subgraph(P,N,graph)
    # stage 9
    graph.merge(s_posts_liked,'User','name')




