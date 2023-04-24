library(rpart)

#' split the data once. Among all N features, randomly sample sqrt(N) features, and choose one as the spliting feature. Selected is based on decision tree. 
#' @param data.train.sub: A data frame with column "y" containing class labels and the remaining columns containing features.
#' @param minsplit: mininum number of genes in a node in order for a split to be attempted. Default value is 2.
#' @param minbucket: mininum number of genes in a node after splitting. Default value is 1.
#' @return A decision tree with only one split.
my.branch = function(data.train.sub, minsplit=2, minbucket=1, verbose=F) {
	features = colnames(data.train.sub); features = features[which(!features %in% c('y', 'idx'))]
	n.features = length(features)
	features.sub = sample(features, ceiling(sqrt(n.features)))
#	features.sub = features
	branch = rpart(y ~ ., data=data.train.sub[, c('y', features.sub)], method='class', control=rpart.control(maxdepth=1, maxcompete=0, maxsurrogate=0, minsplit=2, minbucket=minbucket), model=T, cp=0.01, xval=0)
	if(verbose) cat(features.sub, '... single-level branch generated\n');
	return(branch)
}

#' Build a decision tree with random subsamples of genes and features selected from random subsamples of features. 
#' @param predictors: A data frame with genes in rows and features in columns.
#' @param classes: A vector of class labels.
#' @param random.pct: Perentage to random sample genes.
#' @param minsplit: mininum number of genes in a node in order for a split to be attempted. Default value is 20.
#' @param minbucket: mininum number of genes in a node after splitting. Default value is round(minsplit/3).
#' @return A decision trees.
my.tree = function(predictors, classes, random.pct=0.8, minsplit=20, minbucket=round(minsplit/3), seed=NULL, verbose=F) {
	if(!is.null(seed)) set.seed(seed)
	N = length(classes)
	data.all = data.frame(y=classes, predictors, idx=1:N)
	row.sub = sample(1:N, round(N * random.pct), replace=T) 
	data.train = data.all[row.sub, ]
	data.oob = data.all[-row.sub, ]
	tree = tree.sp = tree.fr = c()
	lev = levels(data.train$y)
	if(verbose) par(mfrow=c(1,5))
	
	in.q = list(list(current=1, parent=0, sample.idx=data.train$idx)); 
	iter = 0;
	flag = 1
	while(flag > 0) {					
		flag = sum(unlist(lapply(in.q, function(x) ifelse(length(x$current) == 0, 0, 1))))
		
		current.node = in.q[[1]]$current
		parent.node = in.q[[1]]$parent
		current.sample = in.q[[1]]$sample.idx
		in.q = in.q[-1]
		iter = iter + 1
		
		data.train.sub = data.train[which(data.train$idx %in% current.sample), ]
		mix = length(which(table(data.train.sub$y) > 0))
		if(length(current.sample) <  minsplit | mix == 1) { 
			in.q[[length(in.q)+1]] = list(current=c(), parent=current.node, sample.idx=c())
			in.q[[length(in.q)+1]] = list(current=c(), parent=current.node, sample.idx=c())
			next;
		}
		branch = my.branch(data.train.sub, minsplit=minsplit, minbucket=minbucket, verbose=verbose)
		if(current.node == 1) {
			tree = branch
			tree$splits = cbind(tree$splits, current.node=current.node, parent.node=parent.node)
		}
		tree.sp = tree$splits; 
		tree.fr = tree$frame;
		
		current.sp = branch$splits; 
		current.fr = branch$frame; 
		if(nrow(current.fr) == 1) { 
			in.q[[length(in.q)+1]] = list(current=c(), parent=current.node, sample.idx=c())
			in.q[[length(in.q)+1]] = list(current=c(), parent=current.node, sample.idx=c())
			next;
		}
		left.idx = right.idx = c()
		if(current.sp[1, 'ncat'] == -1) {
			left.condition = paste('which(data.train.sub$', row.names(current.sp), '<', current.sp[1, 'index'], ')')
			left.idx = eval(parse(text=left.condition))
			right.condition = paste('which(data.train.sub$', row.names(current.sp), '>=', current.sp[1, 'index'], ')')
			right.idx = eval(parse(text=right.condition))
		} else if(current.sp[1, 'ncat'] == 1) {
			left.condition = paste('which(data.train.sub$', row.names(current.sp), '>', current.sp[1, 'index'], ')')
			left.idx = eval(parse(text=left.condition))
			right.condition = paste('which(data.train.sub$', row.names(current.sp), '<=', current.sp[1, 'index'], ')')
			right.idx = eval(parse(text=right.condition))
		} else {
			cat('error: split has wrong ncat - ', current.sp[1, 'ncat'], '\n')
		}		
		## check if left.idx and right.idx match <leaf> n
		if(length(left.idx) != current.fr[2, 'n'] & length(right.idx) != current.fr[3, 'n']) {
			cat('error: left.idx and right.idx have wrong counts.\n')
		}
		
		replace.node = current.node
		if(current.node > 1) {
			child.sp = branch$splits; child.sp = cbind(child.sp, current.node=current.node, parent.node=parent.node)
			child.fr = branch$frame
			yval2 = child.fr$yval2
			if(ncol(yval2) != ncol(tree.fr$yval2)) {
				for(e in 1:length(lev)) {
					if(length(which(data.train.sub$y == lev[e])) == 0) {
						yval2 = cbind(yval2[, 1:e], 0, yval2[, (e+1):(e+2)], 0, yval2[, (e+3):ncol(yval2)]);
						# if(e == 1) {
							# yval2 = cbind(yval2[, 1], 0, yval2[, 2:3], 0, yval2[, 4:ncol(yval2)])
						# } else if(e == 2) {
							# yval2 = cbind(yval2[, 1:2], 0, yval2[, 3:4], 0, yval2[, 5:ncol(yval2)])
						# } else if(e == 3) {
							# yval2 = cbind(yval2[, 1:3], 0, yval2[, 4:5], 0, yval2[, 6:ncol(yval2)])
						# }
					}
				}
				child.fr$yval2 = yval2
			}
			temp.sp = c()
			j = max(which(tree.sp[, 'current.node'] == child.sp[, 'parent.node']), which(tree.sp[, 'parent.node'] == child.sp[, 'parent.node']))
			if(j == nrow(tree.sp)) {
				temp.sp = rbind(tree.sp[1:j, ], child.sp)
				rownames(temp.sp) = c(rownames(tree.sp), rownames(child.sp))
			} else {
				temp.sp = rbind(tree.sp[1:j, ], child.sp, tree.sp[(j+1):nrow(tree.sp), ])
				rownames(temp.sp) = c(rownames(tree.sp)[1:j], rownames(child.sp), rownames(tree.sp)[(j+1):nrow(tree.sp)])
			}
			replace.idx = which(rownames(tree.fr) == current.node)
			temp.fr = c()
			if(replace.idx == nrow(tree.fr)) {
				temp.fr = rbind(tree.fr[1:(replace.idx-1), ], child.fr)
				rownames(temp.fr) = c(rownames(tree.fr), (iter*2):(iter*2+1))
			} else {
				temp.fr = rbind(tree.fr[1:(replace.idx-1), ], child.fr, tree.fr[(replace.idx+1):nrow(tree.fr), ])
				rownames(temp.fr) = c(rownames(tree.fr)[1:(replace.idx)], (iter*2):(iter*2+1), rownames(tree.fr)[(replace.idx+1):nrow(tree.fr)])
			}
			tree$splits = temp.sp
			tree$frame = temp.fr
		}

		in.q[[length(in.q)+1]] = list(current=max(as.numeric(rownames(tree$frame)))-1, parent=replace.node, sample.idx=data.train.sub[left.idx, 'idx'])
		in.q[[length(in.q)+1]] = list(current=max(as.numeric(rownames(tree$frame))), parent=replace.node, sample.idx=data.train.sub[right.idx, 'idx'])

		if(verbose) {
			plot(branch); text(branch);
			plot(tree); text(tree);		
			readline('...');
		}

	}  ## end while
	
	if(nrow(tree$frame) >= 3) {
		tree$terms = terms(as.formula(paste('y ~ ', paste(unique(rownames(tree$splits)), collapse='+'))))
		oob.pred = as.data.frame(predict(tree, newdata=data.oob))
		oob.pred$pred = factor(colnames(oob.pred)[apply(oob.pred, 1, which.max)], levels=c('CIG', 'OG', 'TSG'))
		oob.pred$pred.prob = apply(oob.pred[, 1:3], 1, max)
		tree$data.oob = data.frame(data.oob[, c('idx', 'y')], oob.pred)
	} else {
		tree = NULL
	}
	return(tree)
}

#' Implementation of random forest algorithm. 
#' @param predictors: A data frame with genes in rows and features in columns.
#' @param classes: A vector of class labels.
#' @param random.pct: Percentage to random sample genes.
#' @param ntrees: Number of trees in the random forest. Default value is 200.
#' @param minsplit: mininum number of genes in a node in order for a split to be attempted. Default value is 20.
#' @param minbucket: mininum number of genes in a node after splitting. Default value is round(minsplit/3).
#' @return A list of decision trees.
my.rf = function(predictors, classes, random.pct=0.8, ntrees=200, minsplit=20, minbucket=round(minsplit/3), set.seed=F) {
	forest = list()
	for(ti in 1:ntrees) {
		if(ti %% 50 == 0) cat('Building tree ', ti, '...\n'); flush.console()
		tree = NA
		if(set.seed) {
			tree = my.tree(predictors, classes, random.pct, seed=ti*10, minsplit=minsplit, minbucket=minbucket)
		} else {
			tree = my.tree(predictors, classes, random.pct, minsplit=minsplit, minbucket=minbucket)
		}
		forest[[ti]] = tree
	}
	
	return(forest)
}

#' Make predictions using a random forest. 
#' @param forest: a random forest model.
#' @param newdata: A data frame with genes as rows and features as columns. Prediction is made for each observation.
#' @return A data frame with genes in rows and prediction values in columns. Columns contain prediction probabilities for a gene as a PG, OG, or TSG, respectively, predicted class label with the highest probability.
my.predict = function(forest, newdata) {
	all.pred.list = list()
	for(i in 1:length(forest)) {
		this.tree = forest[[i]]
		if(!is.null(this.tree)) {
			this.pred = as.data.frame(predict(this.tree, newdata=newdata))
			this.pred$pred = factor(colnames(this.pred)[apply(this.pred, 1, which.max)], levels=c('CIG', 'OG', 'TSG'))
			this.pred$pred.prob = apply(this.pred[, 1:3], 1, max)
			this.pred$idx = 1:nrow(this.pred)
			this.pred$tree.idx = i
			all.pred.list[[i]] = this.pred
		}
	}
	all.pred = do.call(rbind, all.pred.list)
	agg = aggregate(tree.idx ~ idx + pred, data=all.pred, 'length')
	agg.pred = do.call(rbind, lapply(split(agg, agg$idx), function(x) {
		temp = data.frame(pred=c('CIG', 'OG', 'TSG'), total=sum(x$tree.idx))
		x = merge(temp, x, by='pred', all=T)
		x[is.na(x)] = 0
		x$prob = x$tree.idx / x$total
		y = c(CIG=x[which(x$pred == 'CIG'), 'prob'], OG=x[which(x$pred == 'OG'), 'prob'], TSG=x[which(x$pred == 'TSG'), 'prob'])
	}))
	agg.pred = data.frame(agg.pred)
	colnames(agg.pred) = c('PG', 'OG', 'TSG')
	agg.pred$pred = colnames(agg.pred)[apply(agg.pred, 1, which.max)]
	agg.pred$pred.prob = apply(agg.pred[, 1:3], 1, max)
	
	pred = data.frame(newdata, agg.pred)
	
	return(pred)
}

#' Calculate performance using out-of-bag samples. 
#' @param forest: a random forest model.
#' @return A list containing two elements. The first element is a data frame with rows correponding to out-of-bag sample and columns correponding to correction predictions. The second element is the overall out-of-bag accuracy.
my.perf = function(forest) {
	all.oob = c()
	for(i in 1:length(forest)) {
		if(!is.null(forest[[i]])) {
			oob = forest[[i]]$data.oob
			oob$correct = ifelse(oob$y == oob$pred, 1, 0)
			oob$tree.idx=i
			all.oob = rbind(all.oob, oob)
		}
	}
	agg.1 = aggregate(tree.idx ~ idx, data=all.oob, 'length')
	agg.2 = aggregate(tree.idx ~ idx, data=all.oob[which(all.oob$correct == 1), ], 'length')
	agg.3 = unique(all.oob[, c('idx', 'y')])
	agg = Reduce(function(x, y) merge(x, y, by='idx', all=T), list(agg.1, agg.2, agg.3))
	colnames(agg) = c('idx', 'cnt.tree', 'cnt.correct', 'y')
	agg[is.na(agg)] = 0
	agg$pct.correct = round(agg$cnt.correct / agg$cnt.tree, 3)
	agg$flag.correct = ifelse(agg$pct.correct > 0.5, 1, 0)
	
	accuracy = round(length(which(agg$flag.correct == 1)) / nrow(agg), 3)
	accuracy = c(all=accuracy, unlist(lapply(split(agg, agg$y), function(x) round(length(which(x$flag.correct == 1)) / nrow(x), 3))))
	
	return(list(agg=agg, accuracy=accuracy));
}

#' Transfer learning by progressive pruning trees in a random forest. 
#' @param forest: a random forest model trained in the source domain.
#' @param newdata: A data frame with genes as rows and features as columns. This is the target domain.
#' @return A list of decision tree.
my.prune = function(tree, newdata) {
	cat('transfer learning by progressive pruning ...\n')
	newdata.scaled = scale(newdata, center=T, scale=T)
	sp = data.frame(tree$splits)
	fr = tree$frame
	in.q = rownames(fr)  ## depth-first traversal
	all.path = path.rpart(tree, in.q, print.it=F);
	all.path = lapply(all.path, function(x) x[-1])	
	all.path.length = unlist(lapply(all.path, length))
	all.subpath = lapply(all.path, function(x) x[-length(x)])
	snip = c()
	while(length(in.q) > 0) {
		this = in.q[1]
		in.q = in.q[-1]		
		this.path = all.path[[this]]
		this.subpath = all.subpath[[this]]
		match.idx = 1:nrow(newdata)  ## root
		child.path.id = c('2', '3') ## child to root
		if(length(this.path) > 0) {
			condition = paste(this.path, collapse='&')
			match.idx = which(with(newdata, eval(parse(text=condition))))
			
			aa = names(all.path.length)[which(all.path.length == length(this.path) + 1)]
			bb = unlist(lapply(all.subpath, function(x) length(intersect(this.path, x))))
			bb = names(all.path.length)[which(bb == length(this.path))]
			child.path.id = intersect(aa, bb)
			if(!length(child.path.id) %in% c(0, 2)) {  ## leave has 0 child; internal node has 2 children
				cat('!! error !! - child.path.id wrong: #', length(child.path.id), '\n'); flush.console()
			}
		}
		if(length(child.path.id) == 0) next;
		
		this.samples = newdata.scaled[match.idx, ]

		child.path.1 = all.path[[child.path.id[1]]]
		condition = paste(child.path.1, collapse='&')
		match.idx = which(with(newdata, eval(parse(text=condition))))
		child.samples.1 = newdata.scaled[match.idx, ]
		
		child.path.2 = all.path[[child.path.id[2]]]
		condition = paste(child.path.2, collapse='&')
		match.idx = which(with(newdata, eval(parse(text=condition))))
		child.samples.2 = newdata.scaled[match.idx, ]
		
		this.dist = mean(dist(this.samples, method='euclidean', diag=F, upper=F))
		child.dist.1 = mean(dist(child.samples.1, method='euclidean', diag=F, upper=F))
		child.dist.2 = mean(dist(child.samples.2, method='euclidean', diag=F, upper=F))
		
		if(is.nan(this.dist) | is.nan(child.dist.1) | is.nan(child.dist.2))  { ## any node that contains 0 samples
			snip = c(snip, this)
		} else if(this.dist < child.dist.1 & this.dist < child.dist.2) {
			snip = c(snip, this)
		}				
	}
	
	cat('snipping ', as.numeric(snip), '\n');
	snipped.tree = tree
	if(length(snip) > 0) {
		snipped.tree = snip.rpart(tree, as.numeric(snip))
	}
	# par(mfrow=c(1,2))
	# plot(tree); text(tree)
	if(length(snipped.tree$splits) > 0) {
		# plot(snipped.tree); text(snipped.tree)
	} else {
		snipped.tree = NULL;
	}
	
	return(snipped.tree)
}

#' Calculate within-node distance after splitting. 
#' @param threshold: threshold used to split the node
#' @param raw.data: Data that will be split.
#' @param norm.data: Data that will be split. Z-transformation has been applied.
#' @param sample.idx: row indices of samples in raw.data and norm.data, for which the splitting is optimized.
#' @param test.var: name of the variable, for which threshold will be shifted.
#' @return mean of within-node distance.
f.dist = function(threshold, raw.data, norm.data, sample.idx, test.var) {
	idx.1 = which(raw.data[, test.var] < threshold)
	idx.2 = which(raw.data[, test.var] >= threshold)
	# idx.1 = intersect(sample.idx, idx.1)
	# idx.2 = intersect(sample.idx, idx.2)
	subset.1 = norm.data[idx.1, ]
	subset.2 = norm.data[idx.2, ]	
	dist.1 = mean(dist(subset.1, method='euclidean', diag=F, upper=F))
	dist.2 = mean(dist(subset.2, method='euclidean', diag=F, upper=F))
	d = mean(c(dist.1, dist.2))
	return(d)
}

#' Transfer learning by threshold shifting trees in a random forest. 
#' @param forest: a random forest model trained in the source domain.
#' @param olddata: A data frame with genes as rows and features as columns. This is the source domain.
#' @param newdata: A data frame with genes as rows and features as columns. This is the target domain.
#' @return A list of decision tree.
my.threshold = function(tree, olddata, newdata) {
	cat('transfer learning by threshold shifting ...\n')
	newdata.scaled = scale(newdata, center=T, scale=T)
	newtree = tree
	sp = newtree$splits
	checked = c()
	th.o_vs_n = c()  ## old threshold vs. new threshold
	for(s in 1:nrow(sp)) {
		this.split = sp[s, ]
		this.parent = sp[s, 'parent.node']
		this.current = sp[s, 'current.node']
		test.var = rownames(sp)[s]
		old.threshold = sp[s, 'index']
		old.percentile = length(which(olddata[, test.var] <= old.threshold))/nrow(olddata)
		new.threshold = quantile(newdata[, test.var], old.percentile) 
#		if(old.threshold * new.threshold < 0)  new.threshold = 0 ## different sign
		if(old.threshold * new.threshold < 0)  new.threshold = old.threshold * 0.5 ## different sign
		old.match.idx = 1:nrow(olddata)  ## root
		new.match.idx = 1:nrow(newdata)  ## root		
		if(this.parent > 0) {
			if(!this.parent %in% checked) {
				cat('!! Error !! parent node', this.parent, ' not checked\n')
			}
			this.path = path.rpart(newtree, this.current, print.it=F)[[1]]
			this.path = this.path[-1]
			condition = paste(this.path, collapse='&')
			old.match.idx = which(with(olddata, eval(parse(text=condition))))
			new.match.idx = which(with(newdata, eval(parse(text=condition))))
		}
		rg = range(c(old.threshold*0.9, old.threshold*1.1, new.threshold*0.9, new.threshold*1.1))
		lower.bound = max(rg[1], min(newdata[, test.var]))
		upper.bound = min(rg[2], max(newdata[, test.var]))
		th = optimize(f.dist, interval=c(lower.bound, upper.bound), raw.data=newdata, norm.data=newdata.scaled, sample.idx=new.match.idx, test.var=test.var)$minimum
		sp[s, 'count'] = length(old.match.idx)
		sp[s, 'index'] = th
		newtree$splits = sp
		checked = c(checked, this.current)
		th.o_vs_n = rbind(th.o_vs_n, c(test.var, old.threshold, th))
	}
	newtree$th.o_vs_n = th.o_vs_n
	
	# par(mfrow=c(1,2))
	# plot(tree); text(tree)
	# plot(newtree); text(newtree)
	
	return(newtree)	
}

#' Transfer learning to adapt a random forest. 
#' @param forest: a random forest model trained in the source domain.
#' @param olddata: A data frame with genes as rows and features as columns. This is the source domain.
#' @param newdata: A data frame with genes as rows and features as columns. This is the target domain.
#' @param method: Adaption via "prune" (progressive pruning), "threshold" (threshold shifting), or "both'.
#' @return A list of decision tree.
my.rf.transfer = function(forest, newdata, olddata, method='both') {
	rf.by_prune = list()
	if(method %in% c('both', 'prune')) {
		for(i in 1:length(forest)) {
			cat('transfer tree', i, '\n'); flush.console();
			this.tree = forest[[i]]
			if(!is.null(this.tree)) {
				snipped.tree = my.prune(this.tree, newdata)
				rf.by_prune[[i]] = snipped.tree
			}
		}
	}
	
	rf.by_threshold = list()
	if(method %in% c('both', 'threshold')) {
		for(i in 1:length(forest)) {
			cat('transfer tree', i, '\n'); flush.console();
			this.tree = forest[[i]]
			if(!is.null(this.tree)) {
				snipped.tree = my.threshold(this.tree, newdata=newdata, olddata=olddata)
				rf.by_threshold[[i]] = snipped.tree
			}
		}
	}

	rf.transfer = c(rf.by_prune, rf.by_threshold)
	return(rf.transfer)
}










