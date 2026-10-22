select
    s.gene_id,
    g.acronym as gene_acronym,
    g.name as gene_name,
    avg(abs(s.mean_lateralisation_index)) as mean_abs_li,
    avg(s.mean_lateralisation_index) as signed_mean_li,
    avg(s.donor_li_stddev) as mean_donor_li_stddev,
    count(*) as region_pair_count
from v_lateralisation_summary s
left join dim_gene g
  on g.gene_id = s.gene_id
group by s.gene_id, g.acronym, g.name
order by mean_abs_li desc
limit 100;
