library(DiagrammeR)

example_graph <-
  create_graph() %>%
  add_pa_graph(
    n = 50, m = 1,
    set_seed = 23
  ) %>%
  add_gnp_graph(
    n = 50, p = 1/100,
    set_seed = 23
  ) %>%
  join_node_attrs(df = get_betweenness(.)) %>%
  join_node_attrs(df = get_degree_total(.)) %>%
  colorize_node_attrs(
    node_attr_from = total_degree,
    node_attr_to = fillcolor,
    palette = "Greens",
    alpha = 90
  ) %>%
  rescale_node_attrs(
    node_attr_from = betweenness,
    to_lower_bound = 0.5,
    to_upper_bound = 1.0,
    node_attr_to = height
  ) %>%
  select_nodes_by_id(nodes = get_articulation_points(.)) %>%
  set_node_attrs_ws(node_attr = peripheries, value = 2) %>%
  set_node_attrs_ws(node_attr = penwidth, value = 3) %>%
  clear_selection() %>%
  set_node_attr_to_display(attr = NULL)

render_graph(example_graph, layout = "nicely")

a_graph <-
  create_graph() %>%
  add_node() %>%
  add_node() %>%
  add_edge(from = 1, to = 2)

render_graph(a_graph, layout = "nicely")

b_graph <- a_graph %>% delete_edge(from = 1, to = 2)
render_graph(b_graph, layout = "nicely")

c_graph <- b_graph %>% add_node(from = 1, to = 2)
render_graph(c_graph, layout = "nicely")

d_graph <-
  c_graph %>%
  add_node(
    type = "type_a",
    node_aes = node_aes(
      color = "steelblue",
      fillcolor = "lightblue",
      fontcolor = "gray35"
    ),
    node_data = node_data(
      value = 2.5
    )
  ) %>%
  add_edge(
    from = 1, to = 3,
    rel = "interacted_with",
    edge_aes = edge_aes(
      color = "red",
      arrowhead = "vee",
      tooltip = "Red Arrow"
    ),
    edge_data = edge_data(
      value = 5.2
    )
  )
render_graph(d_graph, layout = "nicely")

e_graph <-
  d_graph %>%
  select_nodes(conditions = value == 2.5) %>%
  set_node_attrs_ws(node_attr = fillcolor, value = "orange") %>%
  clear_selection()
render_graph(e_graph, layout = "nicely")

f_graph <-
  create_graph() %>%
  add_path(n = 3) %>%
  add_cycle(n = 4) %>%
  add_balanced_tree(k = 2, h = 2)
render_graph(f_graph, layout = "nicely")

render_graph(e_graph, layout = "nicely")

# Create the graph object
i_graph_1 <- create_graph()

# It will start off as empty
i_graph_1 %>% is_graph_empty()

i_graph_2 <-
  i_graph_1 %>%
  add_nodes_from_table(
    table = node_list_1,
    label_col = label
  )

i_graph_2 %>% get_node_df()

i_graph_3 <-
  i_graph_2 %>%
  add_edges_from_table(
    table = edge_list_1,
    from_col = from,
    to_col = to,
    from_to_map = id_external
  )

i_graph_3 %>% get_edge_df()

j_graph <- 
  create_graph() %>% 
  add_nodes_from_table(
    table = node_list_2,
    label_col = label,
    type_col = type
  ) %>%
  add_edges_from_table(
    table = edge_list_2,
    from_col = from,
    to_col = to,
    from_to_map = id_external,
    rel_col = rel
  ) %>%
  drop_node_attrs(node_attr = id_external)

render_graph(j_graph, layout = "nicely")

grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']

      # edge definitions with the node IDs
      tab1 -> tab2 -> tab3 -> tab4 -> tab5;
      }

      [1]: 'Questionnaire sent to n=1000 participants'
      [2]: 'Participants responded to questionnaire n=850'
      [3]: 'Participants came to clinic for evaluation n=700'
      [4]: 'Participants eligible for the study n=600'
      [5]: 'Study sample n=600'
      ")

grViz("digraph flowchart {
      graph [layout = dot]
      # node definitions with substituted label text
      node [fontname = Helvetica, shape = rectangle]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
      tab9 [label = '@@9']
      tab10 [label = '@@10']

      {rank = same; tab6 tab7 tab10}
      # edge definitions with the node IDs
      #tab1 -> tab2 -> tab3 -> tab4 -> tab5;
      tab1 -> tab3 ; tab2 -> tab3 ; tab3 -> tab4 ; tab3 -> tab5 ; tab4 -> tab6 ; tab5 -> tab10 ; tab6 -> tab8; tab7 -> tab8;   tab10 -> tab9 ; tab5 -> tab7
      {rankdir='TB'; tab3 -> tab8 -> tab9}

      #tab8 -> tab9;
      }

      [1]: 'CAF RNA-seq'
      [2]: 'TAN RNA-seq'
      [3]: 'NeoFuse'
      [4]: 'Gene fusions CAF'
      [5]: 'Gene fusions TAN'
      [6]: 'Candidate fusion neoantigens CAF'
      [7]: 'Candidate fusion neoantigens TAN'
      [8]: 'Combined fusion neoantigens'
      [9]: 'Fusion neoantigens specific to CAF'
      [10]: 'Filter out'

      ")
