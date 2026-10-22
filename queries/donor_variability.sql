select
    s.region_pair_id,
    rp.left_name,
    rp.right_name,
    s.gene_id,
    g.acronym as gene_acronym,
    s.mean_lateralisation_index,
    s.donor_li_stddev,
    s.donor_count
from v_lateralisation_summary s
join dim_region_pairs rp
  on rp.region_pair_id = s.region_pair_id
left join dim_gene g
  on g.gene_id = s.gene_id
where s.donor_count > 1
order by s.donor_li_stddev desc
limit 200;
