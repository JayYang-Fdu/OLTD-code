function bk = deinterleaving(ck, alpha)
        [~,s] = sort(alpha);
        bk = ck(s);
end