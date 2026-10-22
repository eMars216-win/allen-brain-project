select
    s.region_pair_id,
    rp.left_name,
    rp.right_name,
    avg(abs(s.mean_lateralisation_index)) as mean_abs_li,
    avg(s.donor_li_stddev) as mean_donor_li_stddev,
    count(*) as gene_count
from v_lateralisation_summary s
join dim_region_pairs rp
  on rp.region_pair_id = s.region_pair_id
group by s.region_pair_id, rp.left_name, rp.right_name
order by mean_abs_li desc;
